#include "EvaluateCM.h"

#include "field.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace OpenOceanKraken
{
    namespace
    {
        using Complex = std::complex<double>;
        using ModComplex = std::complex<float>;

        ModComplex modValue(const Complex value)
        {
            return {static_cast<float>(value.real()), static_cast<float>(value.imag())};
        }

        float interpolateBeamPattern(const SrcBmPat &sbp, const float theta, int &segment)
        {
            while (segment < sbp.NSBPPts - 2 && theta > sbp.theta(segment + 1))
            {
                ++segment;
            }
            while (segment > 0 && theta < sbp.theta(segment))
            {
                --segment;
            }
            const float width = static_cast<float>(sbp.theta(segment + 1) - sbp.theta(segment));
            if (width == 0.0f)
            {
                throw std::runtime_error("Source beam pattern angles must be strictly increasing.");
            }
            const float weight = (theta - static_cast<float>(sbp.theta(segment))) / width;
            return (1.0f - weight) * static_cast<float>(sbp.pat(segment)) +
                   weight * static_cast<float>(sbp.pat(segment + 1));
        }

        float modDepth(const double value)
        {
            return static_cast<float>(value);
        }

        float topDepth(const ssp::SSPStructure &ssp)
        {
            const int medium = ssp.FirstAcoustic;
            if (medium >= 0 && ssp.offset.size() > medium && ssp.z.size() > ssp.offset(medium))
            {
                return modDepth(ssp.z(ssp.offset(medium)));
            }
            return 0.0f;
        }

        float bottomDepth(const ssp::SSPStructure &ssp)
        {
            const int medium = ssp.LastAcoustic;
            const int end = ssp.get_media_end(medium);
            if (end >= 0 && ssp.z.size() > end)
            {
                return modDepth(ssp.z(end));
            }
            const float top = topDepth(ssp);
            return top + (ssp.depth.size() > medium ? modDepth(ssp.depth(medium)) : 0.0f);
        }

        ModComplex halfspaceCp(const HSInfo &halfspace)
        {
            const Complex cp = halfspace.cp != Complex(0.0, 0.0)
                                   ? halfspace.cp
                                   : Complex(halfspace.alphaR, halfspace.alphaI);
            return modValue(cp);
        }

        ModComplex pekerisRoot(const Complex value)
        {
            const Complex root = value.real() >= 0.0 ? std::sqrt(value) : Complex(0.0, 1.0) * std::sqrt(-value);
            return modValue(root);
        }

        ModComplex halfspaceGamma(const ModComplex k, const HSInfo &halfspace, const double frequency)
        {
            const ModComplex cp = halfspaceCp(halfspace);
            const float piSingle = 3.1415926f;
            const ModComplex mediumK = ModComplex(2.0f * piSingle * static_cast<float>(frequency), 0.0f) / cp;
            return pekerisRoot(Complex(k) * Complex(k) - Complex(mediumK) * Complex(mediumK));
        }

        std::vector<float> modeGrid(const EigenParams &eigen)
        {
            std::vector<float> grid(static_cast<size_t>(eigen.ModeZ.size()));
            for (Eigen::Index index = 0; index < eigen.ModeZ.size(); ++index)
            {
                grid[static_cast<size_t>(index)] = modDepth(eigen.ModeZ(index));
            }
            return grid;
        }

        ModComplex modeAtIndex(const EigenParams &eigen, const int mode, const int depthIndex)
        {
            return modValue(eigen.PhiMode(mode, depthIndex));
        }

        ModComplex interpolateMode(const EigenParams &eigen, const int mode, const float depth)
        {
            const std::vector<float> grid = modeGrid(eigen);
            if (depth <= grid.front())
            {
                return modeAtIndex(eigen, mode, 0);
            }
            if (depth >= grid.back())
            {
                return modeAtIndex(eigen, mode, static_cast<int>(grid.size()) - 1);
            }
            const auto upper = std::upper_bound(grid.begin(), grid.end(), depth);
            const int right = static_cast<int>(upper - grid.begin());
            const int left = right - 1;
            const float weight = (depth - grid[static_cast<size_t>(left)]) /
                                 (grid[static_cast<size_t>(right)] - grid[static_cast<size_t>(left)]);
            return modeAtIndex(eigen, mode, left) +
                   weight * (modeAtIndex(eigen, mode, right) - modeAtIndex(eigen, mode, left));
        }

        float mediumBottom(const ssp::SSPStructure &ssp, const int medium)
        {
            const int end = ssp.get_media_end(medium);
            if (end >= 0 && ssp.z.size() > end)
            {
                return modDepth(ssp.z(end));
            }
            const int start = ssp.offset.size() > medium ? ssp.offset(medium) : 0;
            const float top = ssp.z.size() > start ? modDepth(ssp.z(start)) : 0.0f;
            return top + (ssp.depth.size() > medium ? modDepth(ssp.depth(medium)) : 0.0f);
        }

        float densityAt(const ssp::SSPStructure &ssp, const float depth)
        {
            if (depth < topDepth(ssp))
            {
                return static_cast<float>(ssp.HSTop.rho);
            }
            if (depth > bottomDepth(ssp))
            {
                return static_cast<float>(ssp.HSBot.rho);
            }
            for (int medium = ssp.FirstAcoustic; medium <= ssp.LastAcoustic; ++medium)
            {
                if (depth <= mediumBottom(ssp, medium))
                {
                    const int location = ssp.offset.size() > medium ? ssp.offset(medium) : 0;
                    if (ssp.rho.size() > location)
                    {
                        return static_cast<float>(ssp.rho(location));
                    }
                }
            }
            return static_cast<float>(ssp.HSBot.rho);
        }

        float integrationWeight(const ssp::SSPStructure &ssp, const std::vector<float> &grid, const size_t index)
        {
            if (grid.size() < 2)
            {
                throw std::runtime_error("EvaluateCM coupling grid requires at least two depths.");
            }
            if (index == 0)
            {
                const float midpoint = 0.5f * (grid[0] + grid[1]);
                return 0.5f * (grid[1] - grid[0]) / densityAt(ssp, midpoint);
            }
            if (index + 1 == grid.size())
            {
                const float midpoint = 0.5f * (grid[index - 1] + grid[index]);
                return 0.5f * (grid[index] - grid[index - 1]) / densityAt(ssp, midpoint);
            }
            const float leftMidpoint = 0.5f * (grid[index - 1] + grid[index]);
            const float rightMidpoint = 0.5f * (grid[index] + grid[index + 1]);
            return 0.5f * (grid[index] - grid[index - 1]) / densityAt(ssp, leftMidpoint) +
                   0.5f * (grid[index + 1] - grid[index]) / densityAt(ssp, rightMidpoint);
        }

        ModComplex interpolateSamples(const std::vector<float> &grid,
                                      const std::vector<ModComplex> &samples,
                                      const float depth)
        {
            if (depth <= grid.front())
            {
                return samples.front();
            }
            if (depth >= grid.back())
            {
                return samples.back();
            }
            const auto upper = std::upper_bound(grid.begin(), grid.end(), depth);
            const size_t right = static_cast<size_t>(upper - grid.begin());
            const size_t left = right - 1;
            const float weight = (depth - grid[left]) / (grid[right] - grid[left]);
            return samples[left] + weight * (samples[right] - samples[left]);
        }

        ModComplex calculateTail(const float depth,
                                 const std::vector<ModComplex> &leftBoundaryPressure,
                                 const std::vector<ModComplex> &leftGamma,
                                 const float leftBoundaryDepth,
                                 const ModComplex rightMode,
                                 const ModComplex rightGamma,
                                 const float rightBoundaryDepth)
        {
            const ModComplex rightAtDepth = rightMode * std::exp(-rightGamma * (depth - rightBoundaryDepth));
            ModComplex sum(0.0f, 0.0f);
            for (size_t mode = 0; mode < leftBoundaryPressure.size(); ++mode)
            {
                ModComplex numerator = leftBoundaryPressure[mode];
                if (depth != leftBoundaryDepth)
                {
                    numerator *= std::exp(-leftGamma[mode] * (depth - leftBoundaryDepth));
                }
                sum += numerator / (leftGamma[mode] + rightGamma);
            }
            return rightAtDepth * sum;
        }

        void validateInputs(const EigenParams *profiles, const int profileCount,
                            const OOK_parameters &params, const int isz)
        {
            if (profiles == nullptr || params.NProf < 2 || profileCount != params.NProf)
            {
                throw std::invalid_argument("EvaluateCM requires NProf matching at least two eigen profiles.");
            }
            if (params.RProf.size() != params.NProf || params.SSP.size() != static_cast<size_t>(params.NProf) ||
                std::abs(params.RProf(0)) > 1.0e-9)
            {
                throw std::invalid_argument("EvaluateCM requires one SSP and one profile range per profile, starting at zero.");
            }
            for (int profile = 1; profile < params.NProf; ++profile)
            {
                if (!std::isfinite(params.RProf(profile)) || params.RProf(profile) <= params.RProf(profile - 1))
                {
                    throw std::invalid_argument("EvaluateCM profile ranges must be finite and strictly increasing.");
                }
            }
            for (int range = 0; range < params.Pos.NRr; ++range)
            {
                if (!std::isfinite(params.Pos.Rr(range)) || params.Pos.Rr(range) < 0.0 ||
                    (range > 0 && params.Pos.Rr(range) < params.Pos.Rr(range - 1)))
                {
                    throw std::invalid_argument("EvaluateCM receiver ranges must be finite and nondecreasing.");
                }
            }
            if (params.coherenceType != CoherenceType::Coherent)
            {
                throw std::logic_error("EvaluateCM incoherent propagation is not implemented.");
            }
            if (isz < 0 || isz >= params.Pos.NSz)
            {
                throw std::out_of_range("EvaluateCM source index is out of range.");
            }
            for (int profile = 0; profile < profileCount; ++profile)
            {
                const auto &eigen = profiles[profile];
                if (eigen.M <= 0 || eigen.k.size() < eigen.M || eigen.PsiS.rows() < eigen.M ||
                    eigen.PsiS.cols() <= isz || eigen.PsiR.rows() < eigen.M ||
                    eigen.PsiR.cols() < params.Pos.NRz || eigen.ModeZ.size() < 2 ||
                    eigen.PhiMode.rows() < eigen.M || eigen.PhiMode.cols() != eigen.ModeZ.size())
                {
                    throw std::runtime_error("EvaluateCM received an incomplete eigen profile.");
                }
            }
        }

        void advancePhase(std::vector<ModComplex> &amplitudes, const EigenParams &eigen, const double distance)
        {
            for (size_t mode = 0; mode < amplitudes.size(); ++mode)
            {
                const ModComplex k = modValue(eigen.k(static_cast<int>(mode)));
                const Complex phase = std::exp(Complex(0.0, -1.0) * Complex(k) * distance);
                amplitudes[mode] = modValue(Complex(amplitudes[mode]) * phase);
            }
        }

        std::vector<ModComplex> projectToNewProfile(const EigenParams &leftEigen,
                                                     const ssp::SSPStructure &leftSsp,
                                                     const EigenParams &rightEigen,
                                                     const ssp::SSPStructure &rightSsp,
                                                     const std::vector<ModComplex> &leftAmplitudes,
                                                     const double frequency)
        {
            const std::vector<float> leftGrid = modeGrid(leftEigen);
            std::vector<float> couplingGrid = modeGrid(rightEigen);
            for (const float depth : leftGrid)
            {
                if (depth > couplingGrid.back())
                {
                    couplingGrid.push_back(depth);
                }
            }

            const float leftTop = topDepth(leftSsp);
            const float leftBottom = bottomDepth(leftSsp);
            const float rightTop = topDepth(rightSsp);
            const float rightBottom = bottomDepth(rightSsp);
            const bool leftTopOpen = leftSsp.HSTop.BC == BC_Mode::MODE_A_Half_space;
            const bool leftBottomOpen = leftSsp.HSBot.BC == BC_Mode::MODE_A_Half_space;
            const bool rightTopOpen = rightSsp.HSTop.BC == BC_Mode::MODE_A_Half_space;
            const bool rightBottomOpen = rightSsp.HSBot.BC == BC_Mode::MODE_A_Half_space;

            std::vector<ModComplex> leftPressure(leftGrid.size(), ModComplex(0.0f, 0.0f));
            std::vector<ModComplex> leftTopPressure(leftAmplitudes.size());
            std::vector<ModComplex> leftBottomPressure(leftAmplitudes.size());
            std::vector<ModComplex> leftTopGamma(leftAmplitudes.size(), ModComplex(0.0f, 0.0f));
            std::vector<ModComplex> leftBottomGamma(leftAmplitudes.size(), ModComplex(0.0f, 0.0f));
            for (size_t mode = 0; mode < leftAmplitudes.size(); ++mode)
            {
                for (size_t depth = 0; depth < leftGrid.size(); ++depth)
                {
                    leftPressure[depth] += leftAmplitudes[mode] *
                                           modeAtIndex(leftEigen, static_cast<int>(mode), static_cast<int>(depth));
                }
                leftTopPressure[mode] = leftAmplitudes[mode] * modeAtIndex(leftEigen, static_cast<int>(mode), 0);
                leftBottomPressure[mode] = leftAmplitudes[mode] *
                    modeAtIndex(leftEigen, static_cast<int>(mode), static_cast<int>(leftGrid.size()) - 1);
                const ModComplex k = modValue(leftEigen.k(static_cast<int>(mode)));
                if (leftTopOpen)
                {
                    leftTopGamma[mode] = halfspaceGamma(k, leftSsp.HSTop, frequency);
                }
                if (leftBottomOpen)
                {
                    leftBottomGamma[mode] = halfspaceGamma(k, leftSsp.HSBot, frequency);
                }
            }

            std::vector<ModComplex> weightedPressure(couplingGrid.size(), ModComplex(0.0f, 0.0f));
            for (size_t depthIndex = 0; depthIndex < couplingGrid.size(); ++depthIndex)
            {
                const float depth = couplingGrid[depthIndex];
                ModComplex pressure(0.0f, 0.0f);
                if (depth > leftBottom)
                {
                    if (leftBottomOpen)
                    {
                        for (size_t mode = 0; mode < leftAmplitudes.size(); ++mode)
                        {
                            pressure += leftBottomPressure[mode] *
                                        std::exp(-leftBottomGamma[mode] * (depth - leftBottom));
                        }
                    }
                }
                else if (depth < leftTop)
                {
                    if (leftTopOpen)
                    {
                        for (size_t mode = 0; mode < leftAmplitudes.size(); ++mode)
                        {
                            pressure += leftTopPressure[mode] *
                                        std::exp(-leftTopGamma[mode] * (leftTop - depth));
                        }
                    }
                }
                else
                {
                    pressure = interpolateSamples(leftGrid, leftPressure, depth);
                }
                weightedPressure[depthIndex] = integrationWeight(rightSsp, couplingGrid, depthIndex) * pressure;
            }

            std::vector<ModComplex> rightAmplitudes(static_cast<size_t>(rightEigen.M), ModComplex(0.0f, 0.0f));
            const float rightTopDensity = static_cast<float>(rightSsp.HSTop.rho);
            const float rightBottomDensity = static_cast<float>(rightSsp.HSBot.rho);
            for (int mode = 0; mode < rightEigen.M; ++mode)
            {
                const ModComplex k = modValue(rightEigen.k(mode));
                const ModComplex topMode = modeAtIndex(rightEigen, mode, 0);
                const ModComplex bottomMode = modeAtIndex(rightEigen, mode, rightEigen.PhiMode.cols() - 1);
                const ModComplex topGamma = rightTopOpen ? halfspaceGamma(k, rightSsp.HSTop, frequency) : ModComplex(0.0f, 0.0f);
                const ModComplex bottomGamma = rightBottomOpen ? halfspaceGamma(k, rightSsp.HSBot, frequency) : ModComplex(0.0f, 0.0f);

                ModComplex projection(0.0f, 0.0f);
                for (size_t depthIndex = 0; depthIndex < couplingGrid.size(); ++depthIndex)
                {
                    const float depth = couplingGrid[depthIndex];
                    ModComplex rightMode(0.0f, 0.0f);
                    if (depth > rightBottom)
                    {
                        if (rightBottomOpen)
                        {
                            rightMode = bottomMode * std::exp(-bottomGamma * (depth - rightBottom));
                        }
                    }
                    else if (depth < rightTop)
                    {
                        if (rightTopOpen)
                        {
                            rightMode = topMode * std::exp(-topGamma * (rightTop - depth));
                        }
                    }
                    else
                    {
                        rightMode = interpolateMode(rightEigen, mode, depth);
                    }
                    projection += weightedPressure[depthIndex] * rightMode;
                }

                if (rightTopOpen)
                {
                    projection += calculateTail(couplingGrid.front(), leftTopPressure, leftTopGamma, leftTop,
                                                topMode / rightTopDensity, topGamma, rightTop);
                }
                if (rightBottomOpen)
                {
                    projection += calculateTail(couplingGrid.back(), leftBottomPressure, leftBottomGamma, leftBottom,
                                                bottomMode / rightBottomDensity, bottomGamma, rightBottom);
                }
                rightAmplitudes[static_cast<size_t>(mode)] = projection;
            }
            return rightAmplitudes;
        }
    }

    void EvaluateCM(const EigenParams *profiles,
                    const int profileCount,
                    const OOK_parameters &params,
                    const int isz,
                    std::complex<float> *uAllSources)
    {
        validateInputs(profiles, profileCount, params, isz);
        if (uAllSources == nullptr)
        {
            throw std::invalid_argument("EvaluateCM pressure output is null.");
        }

        int profile = 0;
        const int sourceModes = std::min(params.MLimit, profiles[0].M);
        if (sourceModes <= 0)
        {
            throw std::runtime_error("EvaluateCM has no propagating source modes after MLimit.");
        }

        constexpr float piSingle = 3.1415926f;
        const ModComplex imaginary(0.0f, 1.0f);
        ModComplex factor = std::sqrt(2.0f * piSingle) * std::exp(imaginary * (piSingle / 4.0f));
        if (params.SourceType == Source_Mode::MODE_R_Point)
        {
            factor *= imaginary;
        }

        std::vector<ModComplex> amplitudes(static_cast<size_t>(sourceModes));
        int beamSegment = 0;
        for (int mode = 0; mode < sourceModes; ++mode)
        {
            ModComplex source = modValue(profiles[0].PsiS(mode, isz));
            const ModComplex k = modValue(profiles[0].k(mode));
            if (params.SBP.isSet && isz == 0 && params.SBP.NSBPPts >= 2)
            {
                const float omega = 2.0f * piSingle * static_cast<float>(params.freqinfo.freq);
                const float kz2 = std::max(0.0f, std::real(omega * omega / (1500.0f * 1500.0f) - k * k));
                const float theta = (180.0f / piSingle) * std::atan(std::sqrt(kz2) / std::real(k));
                source *= interpolateBeamPattern(params.SBP, theta, beamSegment);
            }
            amplitudes[static_cast<size_t>(mode)] = params.SourceType == Source_Mode::MODE_X_Line
                ? factor * source / k
                : factor * source / std::sqrt(k);
        }

        std::vector<double> interfaces(static_cast<size_t>(params.NProf), 0.0);
        for (int index = 1; index < params.NProf; ++index)
        {
            interfaces[static_cast<size_t>(index)] = 0.5 * (params.RProf(index - 1) + params.RProf(index));
        }

        for (int ir = 0; ir < params.Pos.NRr; ++ir)
        {
            const double receiverRange = params.Pos.Rr(ir);
            const double previousRange = ir > 0 ? params.Pos.Rr(ir - 1) : 0.0;
            if (profile + 1 < params.NProf && receiverRange > interfaces[static_cast<size_t>(profile + 1)])
            {
                ++profile;
                advancePhase(amplitudes, profiles[profile - 1],
                             interfaces[static_cast<size_t>(profile)] - previousRange);
                amplitudes = projectToNewProfile(profiles[profile - 1], params.SSP[static_cast<size_t>(profile - 1)],
                                                 profiles[profile], params.SSP[static_cast<size_t>(profile)],
                                                 amplitudes, params.freqinfo.freq);

                while (profile + 1 < params.NProf && receiverRange > interfaces[static_cast<size_t>(profile + 1)])
                {
                    ++profile;
                    advancePhase(amplitudes, profiles[profile - 1],
                                 interfaces[static_cast<size_t>(profile)] - interfaces[static_cast<size_t>(profile - 1)]);
                    amplitudes = projectToNewProfile(profiles[profile - 1], params.SSP[static_cast<size_t>(profile - 1)],
                                                     profiles[profile], params.SSP[static_cast<size_t>(profile)],
                                                     amplitudes, params.freqinfo.freq);
                }
                advancePhase(amplitudes, profiles[profile],
                             receiverRange - interfaces[static_cast<size_t>(profile)]);
            }
            else
            {
                advancePhase(amplitudes, profiles[profile], receiverRange - previousRange);
            }

            for (int iz = 0; iz < params.Pos.NRz; ++iz)
            {
                ModComplex pressure(0.0f, 0.0f);
                for (size_t mode = 0; mode < amplitudes.size(); ++mode)
                {
                    pressure += amplitudes[mode] * modValue(profiles[profile].PsiR(static_cast<int>(mode), iz));
                }
                if (params.SourceType == Source_Mode::MODE_R_Point && receiverRange != 0.0)
                {
                    pressure /= std::sqrt(static_cast<float>(receiverRange));
                }
                uAllSources[GetFieldAddr(isz, iz, ir, &params.Pos)] = pressure;
            }
        }
    }
}
