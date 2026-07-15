#include "EvaluateAD.h"

#include "field.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace OpenOceanKraken
{
    namespace
    {
        using Complex = std::complex<double>;
        using ModComplex = std::complex<float>;

        Complex roundMod(const Complex value)
        {
            const ModComplex rounded(static_cast<float>(value.real()), static_cast<float>(value.imag()));
            return {static_cast<double>(rounded.real()), static_cast<double>(rounded.imag())};
        }

        Complex interpolateMod(const Complex left, const Complex right, const double weight)
        {
            return roundMod(roundMod(left) + weight * (roundMod(right) - roundMod(left)));
        }

        double interpolateBeamPattern(const SrcBmPat &sbp, const double theta, int &segment)
        {
            const int count = sbp.NSBPPts;
            while (segment < count - 2 && theta > sbp.theta(segment + 1))
            {
                ++segment;
            }
            while (segment > 0 && theta < sbp.theta(segment))
            {
                --segment;
            }
            const double width = sbp.theta(segment + 1) - sbp.theta(segment);
            if (width == 0.0)
            {
                throw std::runtime_error("Source beam pattern angles must be strictly increasing.");
            }
            const double weight = (theta - sbp.theta(segment)) / width;
            return (1.0 - weight) * sbp.pat(segment) + weight * sbp.pat(segment + 1);
        }

        void validateInputs(const EigenParams *profiles, const int profileCount,
                            const OOK_parameters &params, const int isz)
        {
            if (profiles == nullptr || params.NProf < 2 || params.NProf != profileCount)
            {
                throw std::invalid_argument("EvaluateAD requires NProf matching at least two eigen profiles.");
            }
            if (params.RProf.size() != params.NProf || std::abs(params.RProf(0)) > 1.0e-9)
            {
                throw std::invalid_argument("EvaluateAD requires a complete RProf vector starting at zero.");
            }
            if (isz < 0 || isz >= params.Pos.NSz)
            {
                throw std::out_of_range("EvaluateAD source index is out of range.");
            }
            for (int profile = 0; profile < params.NProf; ++profile)
            {
                const auto &eigen = profiles[profile];
                if (eigen.M <= 0 || eigen.k.size() < eigen.M || eigen.PsiS.rows() < eigen.M ||
                    eigen.PsiS.cols() <= isz || eigen.PsiR.rows() < eigen.M || eigen.PsiR.cols() < params.Pos.NRz)
                {
                    throw std::runtime_error("EvaluateAD received an incomplete eigen profile.");
                }
                if (profile > 0 && params.RProf(profile) <= params.RProf(profile - 1))
                {
                    throw std::invalid_argument("EvaluateAD profile ranges must be strictly increasing.");
                }
            }
            for (int range = 0; range < params.Pos.NRr; ++range)
            {
                if (!std::isfinite(params.Pos.Rr(range)) || params.Pos.Rr(range) < 0.0 ||
                    (range > 0 && params.Pos.Rr(range) < params.Pos.Rr(range - 1)))
                {
                    throw std::invalid_argument("EvaluateAD receiver ranges must be finite and nondecreasing.");
                }
            }
        }
    }

    void EvaluateAD(const EigenParams *profiles,
                    const int profileCount,
                    const OOK_parameters &params,
                    const int isz,
                    std::complex<float> *uAllSources)
    {
        validateInputs(profiles, profileCount, params, isz);
        if (uAllSources == nullptr)
        {
            throw std::invalid_argument("EvaluateAD pressure output is null.");
        }

        int modeCount = std::min({params.MLimit, profiles[0].M, profiles[1].M});
        if (modeCount <= 0)
        {
            throw std::runtime_error("EvaluateAD has no propagating modes after MLimit.");
        }

        constexpr float piSingle = 3.1415926f;
        const ModComplex imaginary(0.0f, 1.0f);
        const ModComplex factorSingle = imaginary * std::sqrt(2.0f * piSingle) *
                                        std::exp(imaginary * (piSingle / 4.0f));

        std::vector<Complex> constants(static_cast<size_t>(modeCount));
        for (int mode = 0; mode < modeCount; ++mode)
        {
            ModComplex source(static_cast<float>(profiles[0].PsiS(mode, isz).real()),
                              static_cast<float>(profiles[0].PsiS(mode, isz).imag()));
            constants[static_cast<size_t>(mode)] = roundMod(factorSingle * source);
        }

        if (params.SBP.isSet && isz == 0 && params.SBP.NSBPPts >= 2)
        {
            int segment = 0;
            const double omega = 2.0 * 3.14159265358979323846 * params.freqinfo.freq;
            for (int mode = 0; mode < modeCount; ++mode)
            {
                const Complex wavenumber = roundMod(profiles[0].k(mode));
                const double kz2 = std::max(0.0, std::real(omega * omega / (1500.0 * 1500.0) - wavenumber * wavenumber));
                const double theta = (180.0 / 3.14159265358979323846) *
                                     std::atan(std::sqrt(kz2) / std::real(wavenumber));
                constants[static_cast<size_t>(mode)] *= interpolateBeamPattern(params.SBP, theta, segment);
            }
        }

        std::vector<Complex> sumK(static_cast<size_t>(modeCount), Complex(0.0, 0.0));
        std::vector<Complex> sumKInv(static_cast<size_t>(modeCount), Complex(0.0, 0.0));
        int profile = 0;

        auto profileEnd = [&params](const int index)
        {
            return index + 1 < params.NProf ? params.RProf(index + 1) : 1.0e23;
        };
        auto interpolatedK = [profiles, &params](const int index, const int mode, const double weight)
        {
            const int right = std::min(index + 1, params.NProf - 1);
            return interpolateMod(profiles[index].k(mode), profiles[right].k(mode), weight);
        };

        for (int ir = 0; ir < params.Pos.NRr; ++ir)
        {
            const double receiverRange = params.Pos.Rr(ir);
            while (profile + 1 < params.NProf && receiverRange > profileEnd(profile))
            {
                const double left = ir > 0 ? std::max(params.Pos.Rr(ir - 1), params.RProf(profile)) : params.RProf(profile);
                const double right = profileEnd(profile);
                const double midpoint = 0.5 * (right + left);
                const double midpointWeight = (midpoint - params.RProf(profile)) / (right - params.RProf(profile));
                for (int mode = 0; mode < modeCount; ++mode)
                {
                    const Complex kMid = interpolatedK(profile, mode, midpointWeight);
                    sumK[static_cast<size_t>(mode)] += kMid * (right - left);
                    sumKInv[static_cast<size_t>(mode)] += (right - left) / kMid;
                }
                ++profile;
                if (profile + 1 < params.NProf)
                {
                    modeCount = std::min(modeCount, profiles[profile + 1].M);
                }
            }

            const double left = ir > 0 ? std::max(params.Pos.Rr(ir - 1), params.RProf(profile)) : params.RProf(profile);
            const double right = profileEnd(profile);
            const double midpoint = 0.5 * (receiverRange + left);
            const double weight = (receiverRange - params.RProf(profile)) / (right - params.RProf(profile));
            const double midpointWeight = (midpoint - params.RProf(profile)) / (right - params.RProf(profile));

            std::vector<Complex> hank(static_cast<size_t>(modeCount));
            std::vector<Complex> kAtRange(static_cast<size_t>(modeCount));
            for (int mode = 0; mode < modeCount; ++mode)
            {
                const Complex kMid = interpolatedK(profile, mode, midpointWeight);
                const Complex kInt = interpolatedK(profile, mode, weight);
                kAtRange[static_cast<size_t>(mode)] = kInt;
                sumK[static_cast<size_t>(mode)] += kMid * (receiverRange - left);
                sumKInv[static_cast<size_t>(mode)] += (receiverRange - left) / kMid;

                Complex propagation;
                if (params.coherenceType == CoherenceType::Coherent)
                {
                    propagation = roundMod(std::exp(Complex(0.0, -1.0) * sumK[static_cast<size_t>(mode)]));
                }
                else
                {
                    propagation = std::exp(std::real(Complex(0.0, -1.0) * sumK[static_cast<size_t>(mode)]));
                }
                hank[static_cast<size_t>(mode)] = constants[static_cast<size_t>(mode)] * propagation;

                if (params.SourceType == Source_Mode::MODE_R_Point)
                {
                    hank[static_cast<size_t>(mode)] = receiverRange == 0.0
                        ? Complex(0.0, 0.0)
                        : hank[static_cast<size_t>(mode)] / std::sqrt(kInt * receiverRange);
                }
                else
                {
                    hank[static_cast<size_t>(mode)] /= kInt;
                }
            }

            const int rightProfile = std::min(profile + 1, params.NProf - 1);
            for (int iz = 0; iz < params.Pos.NRz; ++iz)
            {
                Complex pressure(0.0, 0.0);
                double incoherentEnergy = 0.0;
                for (int mode = 0; mode < modeCount; ++mode)
                {
                    const Complex phi = interpolateMod(
                        profiles[profile].PsiR(mode, iz),
                        profiles[rightProfile].PsiR(mode, iz), weight);
                    const Complex term = phi * hank[static_cast<size_t>(mode)];
                    if (params.coherenceType == CoherenceType::Coherent)
                    {
                        pressure += term;
                    }
                    else
                    {
                        incoherentEnergy += std::norm(term);
                    }
                }
                if (params.coherenceType == CoherenceType::Incoherent)
                {
                    pressure = {std::sqrt(incoherentEnergy), 0.0};
                }
                uAllSources[GetFieldAddr(isz, iz, ir, &params.Pos)] = ModComplex(
                    static_cast<float>(pressure.real()), static_cast<float>(pressure.imag()));
            }
        }
    }
}
