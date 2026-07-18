#include "algorithm/FieldSolver.h"

#include "algorithm/ComplexNumerics.h"

#include <algorithm>
#include <cmath>
#include <future>
#include <stdexcept>
#include <vector>

namespace OpenOceanKrakenc
{
namespace
{
using ModComplex = std::complex<float>;

ModComplex modValue(std::complex<double> value)
{
    return {static_cast<float>(value.real()), static_cast<float>(value.imag())};
}

float topDepth(const ModeProfileData &profile)
{
    return static_cast<float>(profile.top.depth);
}

float bottomDepth(const ModeProfileData &profile)
{
    return static_cast<float>(profile.bottom.depth);
}

float densityAt(const ModeProfileData &profile, float depth)
{
    if (depth < topDepth(profile))
    {
        return static_cast<float>(profile.top.rho);
    }
    if (depth > bottomDepth(profile))
    {
        return static_cast<float>(profile.bottom.rho);
    }
    for (std::size_t medium = 0; medium < profile.mediumDepths.size(); ++medium)
    {
        const float nextDepth = medium + 1 < profile.mediumDepths.size()
            ? static_cast<float>(profile.mediumDepths[medium + 1])
            : bottomDepth(profile);
        if (depth <= nextDepth)
        {
            return static_cast<float>(profile.mediumDensities[medium]);
        }
    }
    return static_cast<float>(profile.bottom.rho);
}

float integrationWeight(const ModeProfileData &profile,
                        const std::vector<float> &grid,
                        std::size_t index)
{
    if (grid.size() < 2)
    {
        throw std::runtime_error("coupled-mode integration grid is too short");
    }
    if (index == 0)
    {
        const float midpoint = 0.5f * (grid[0] + grid[1]);
        return 0.5f * (grid[1] - grid[0]) / densityAt(profile, midpoint);
    }
    if (index + 1 == grid.size())
    {
        const float midpoint = 0.5f * (grid[index - 1] + grid[index]);
        return 0.5f * (grid[index] - grid[index - 1]) /
               densityAt(profile, midpoint);
    }
    const float leftMidpoint = 0.5f * (grid[index - 1] + grid[index]);
    const float rightMidpoint = 0.5f * (grid[index] + grid[index + 1]);
    return 0.5f * (grid[index] - grid[index - 1]) /
               densityAt(profile, leftMidpoint) +
           0.5f * (grid[index + 1] - grid[index]) /
               densityAt(profile, rightMidpoint);
}

ModComplex halfspaceGamma(ModComplex k,
                          const ModeBoundaryData &boundary,
                          double frequency)
{
    constexpr float piSingle = 3.1415926f;
    const ModComplex cp = modValue(boundary.cp);
    const ModComplex mediumK =
        ModComplex(2.0f * piSingle * static_cast<float>(frequency), 0.0f) / cp;
    return modValue(pekerisRoot(std::complex<double>(k) * std::complex<double>(k) -
                                std::complex<double>(mediumK) *
                                    std::complex<double>(mediumK)));
}

ModComplex modeAtGrid(const ModeProfileData &profile,
                      std::size_t mode,
                      std::size_t depth)
{
    return modValue(profile.modes[mode][depth]);
}

ModComplex interpolateSamples(const std::vector<float> &grid,
                              const std::vector<ModComplex> &samples,
                              float depth)
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
    const std::size_t right = static_cast<std::size_t>(upper - grid.begin());
    const std::size_t left = right - 1;
    const float weight = (depth - grid[left]) / (grid[right] - grid[left]);
    return samples[left] + weight * (samples[right] - samples[left]);
}

ModComplex modeAtDepth(const ModeProfileData &profile,
                       std::size_t mode,
                       float depth,
                       double frequency)
{
    if (depth < topDepth(profile))
    {
        if (profile.top.type != 'A')
        {
            return modeAtGrid(profile, mode, 0);
        }
        const ModComplex gamma = halfspaceGamma(
            modValue(profile.wavenumbers[mode]), profile.top, frequency);
        return modeAtGrid(profile, mode, 0) *
               std::exp(-gamma * (topDepth(profile) - depth));
    }
    if (depth > bottomDepth(profile))
    {
        if (profile.bottom.type != 'A')
        {
            return modeAtGrid(profile, mode, profile.depths.size() - 1);
        }
        const ModComplex gamma = halfspaceGamma(
            modValue(profile.wavenumbers[mode]), profile.bottom, frequency);
        return modeAtGrid(profile, mode, profile.depths.size() - 1) *
               std::exp(-gamma * (depth - bottomDepth(profile)));
    }
    if (depth <= profile.depths.front())
    {
        return modeAtGrid(profile, mode, 0);
    }
    if (depth >= profile.depths.back())
    {
        return modeAtGrid(profile, mode, profile.depths.size() - 1);
    }
    const auto upper = std::upper_bound(profile.depths.begin(), profile.depths.end(),
                                        static_cast<double>(depth));
    const std::size_t right = static_cast<std::size_t>(upper - profile.depths.begin());
    const std::size_t left = right - 1;
    const float weight = (depth - static_cast<float>(profile.depths[left])) /
        static_cast<float>(profile.depths[right] - profile.depths[left]);
    return modeAtGrid(profile, mode, left) + weight *
        (modeAtGrid(profile, mode, right) - modeAtGrid(profile, mode, left));
}

ModComplex modeDerivativeAtDepth(const ModeProfileData &profile,
                                 std::size_t mode,
                                 float depth,
                                 double frequency)
{
    if (depth < topDepth(profile) && profile.top.type == 'A')
    {
        return halfspaceGamma(modValue(profile.wavenumbers[mode]),
                              profile.top, frequency) *
               modeAtDepth(profile, mode, depth, frequency);
    }
    if (depth > bottomDepth(profile) && profile.bottom.type == 'A')
    {
        return -halfspaceGamma(modValue(profile.wavenumbers[mode]),
                               profile.bottom, frequency) *
               modeAtDepth(profile, mode, depth, frequency);
    }
    if (profile.depths.size() < 2)
    {
        return {0.0f, 0.0f};
    }
    std::size_t right = 1;
    if (depth >= profile.depths.back())
    {
        right = profile.depths.size() - 1;
    }
    else if (depth > profile.depths.front())
    {
        right = static_cast<std::size_t>(std::upper_bound(
            profile.depths.begin(), profile.depths.end(),
            static_cast<double>(depth)) - profile.depths.begin());
    }
    const std::size_t left = right - 1;
    return (modeAtGrid(profile, mode, right) - modeAtGrid(profile, mode, left)) /
           static_cast<float>(profile.depths[right] - profile.depths[left]);
}

ModComplex calculateTail(float depth,
                         const std::vector<ModComplex> &leftBoundaryPressure,
                         const std::vector<ModComplex> &leftGamma,
                         float leftBoundaryDepth,
                         ModComplex rightMode,
                         ModComplex rightGamma,
                         float rightBoundaryDepth)
{
    const ModComplex rightAtDepth =
        rightMode * std::exp(-rightGamma * (depth - rightBoundaryDepth));
    ModComplex sum(0.0f, 0.0f);
    for (std::size_t mode = 0; mode < leftBoundaryPressure.size(); ++mode)
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

void advancePhase(std::vector<ModComplex> &amplitudes,
                  const ModeProfileData &profile,
                  double distance)
{
    for (std::size_t mode = 0; mode < amplitudes.size(); ++mode)
    {
        const ModComplex k = modValue(profile.wavenumbers[mode]);
        const std::complex<double> phase =
            std::exp(std::complex<double>(0.0, -1.0) *
                     std::complex<double>(k) * distance);
        amplitudes[mode] = modValue(std::complex<double>(amplitudes[mode]) * phase);
    }
}

std::vector<ModComplex> projectToNewProfile(
    const ModeProfileData &left,
    const ModeProfileData &right,
    const std::vector<ModComplex> &leftAmplitudes,
    double frequency)
{
    std::vector<float> leftGrid(left.depths.size());
    std::vector<float> couplingGrid(right.depths.size());
    std::transform(left.depths.begin(), left.depths.end(), leftGrid.begin(),
                   [](double value) { return static_cast<float>(value); });
    std::transform(right.depths.begin(), right.depths.end(), couplingGrid.begin(),
                   [](double value) { return static_cast<float>(value); });
    for (float depth : leftGrid)
    {
        if (depth > couplingGrid.back())
        {
            couplingGrid.push_back(depth);
        }
    }

    const bool leftTopOpen = left.top.type == 'A';
    const bool leftBottomOpen = left.bottom.type == 'A';
    const bool rightTopOpen = right.top.type == 'A';
    const bool rightBottomOpen = right.bottom.type == 'A';
    std::vector<ModComplex> leftPressure(leftGrid.size(), ModComplex(0.0f, 0.0f));
    std::vector<ModComplex> leftTopPressure(leftAmplitudes.size());
    std::vector<ModComplex> leftBottomPressure(leftAmplitudes.size());
    std::vector<ModComplex> leftTopGamma(leftAmplitudes.size(), ModComplex(0.0f, 0.0f));
    std::vector<ModComplex> leftBottomGamma(leftAmplitudes.size(), ModComplex(0.0f, 0.0f));
    for (std::size_t mode = 0; mode < leftAmplitudes.size(); ++mode)
    {
        for (std::size_t depth = 0; depth < leftGrid.size(); ++depth)
        {
            leftPressure[depth] += leftAmplitudes[mode] *
                                   modeAtGrid(left, mode, depth);
        }
        leftTopPressure[mode] = leftAmplitudes[mode] * modeAtGrid(left, mode, 0);
        leftBottomPressure[mode] = leftAmplitudes[mode] *
            modeAtGrid(left, mode, leftGrid.size() - 1);
        const ModComplex k = modValue(left.wavenumbers[mode]);
        if (leftTopOpen)
        {
            leftTopGamma[mode] = halfspaceGamma(k, left.top, frequency);
        }
        if (leftBottomOpen)
        {
            leftBottomGamma[mode] = halfspaceGamma(k, left.bottom, frequency);
        }
    }

    std::vector<ModComplex> weightedPressure(couplingGrid.size(), ModComplex(0.0f, 0.0f));
    for (std::size_t depthIndex = 0; depthIndex < couplingGrid.size(); ++depthIndex)
    {
        const float depth = couplingGrid[depthIndex];
        ModComplex pressure(0.0f, 0.0f);
        if (depth > bottomDepth(left))
        {
            if (leftBottomOpen)
            {
                for (std::size_t mode = 0; mode < leftAmplitudes.size(); ++mode)
                {
                    pressure += leftBottomPressure[mode] * std::exp(
                        -leftBottomGamma[mode] * (depth - bottomDepth(left)));
                }
            }
        }
        else if (depth < topDepth(left))
        {
            if (leftTopOpen)
            {
                for (std::size_t mode = 0; mode < leftAmplitudes.size(); ++mode)
                {
                    pressure += leftTopPressure[mode] * std::exp(
                        -leftTopGamma[mode] * (topDepth(left) - depth));
                }
            }
        }
        else
        {
            pressure = interpolateSamples(leftGrid, leftPressure, depth);
        }
        weightedPressure[depthIndex] =
            integrationWeight(right, couplingGrid, depthIndex) * pressure;
    }

    std::vector<ModComplex> rightAmplitudes(
        right.modes.size(), ModComplex(0.0f, 0.0f));
    for (std::size_t mode = 0; mode < right.modes.size(); ++mode)
    {
        const ModComplex k = modValue(right.wavenumbers[mode]);
        const ModComplex topMode = modeAtGrid(right, mode, 0);
        const ModComplex bottomMode = modeAtGrid(right, mode, right.depths.size() - 1);
        const ModComplex topGamma = rightTopOpen
            ? halfspaceGamma(k, right.top, frequency) : ModComplex(0.0f, 0.0f);
        const ModComplex bottomGamma = rightBottomOpen
            ? halfspaceGamma(k, right.bottom, frequency) : ModComplex(0.0f, 0.0f);
        ModComplex projection(0.0f, 0.0f);
        for (std::size_t depthIndex = 0; depthIndex < couplingGrid.size(); ++depthIndex)
        {
            projection += weightedPressure[depthIndex] *
                          modeAtDepth(right, mode, couplingGrid[depthIndex], frequency);
        }
        if (rightTopOpen)
        {
            projection += calculateTail(
                couplingGrid.front(), leftTopPressure, leftTopGamma, topDepth(left),
                topMode / static_cast<float>(right.top.rho), topGamma, topDepth(right));
        }
        if (rightBottomOpen)
        {
            projection += calculateTail(
                couplingGrid.back(), leftBottomPressure, leftBottomGamma,
                bottomDepth(left), bottomMode / static_cast<float>(right.bottom.rho),
                bottomGamma, bottomDepth(right));
        }
        rightAmplitudes[mode] = projection;
    }
    return rightAmplitudes;
}

template <typename Function>
void parallelFor(std::size_t count, std::size_t requestedThreads,
                 const Function &function)
{
    const std::size_t workers = std::min(
        std::max<std::size_t>(1, requestedThreads), count);
    if (workers == 1)
    {
        for (std::size_t index = 0; index < count; ++index)
        {
            function(index);
        }
        return;
    }
    for (std::size_t batch = 0; batch < count; batch += workers)
    {
        std::vector<std::future<void>> futures;
        const std::size_t end = std::min(count, batch + workers);
        for (std::size_t index = batch; index < end; ++index)
        {
            futures.push_back(std::async(std::launch::async, function, index));
        }
        for (auto &future : futures)
        {
            future.get();
        }
    }
}
}

PressureField evaluateCoupledField(
    const ModeFileData &modes,
    const FieldParameters &parameters)
{
    validateFieldParameters(parameters);
    std::vector<const ModeProfileData *> profiles{&modes};
    for (const ModeProfileData &profile : modes.additionalProfiles)
    {
        profiles.push_back(&profile);
    }
    if (profiles.size() < 2 || parameters.profileRangesKm.size() != profiles.size() ||
        parameters.profileRangesKm.front() != 0.0 || !parameters.coherent)
    {
        throw std::invalid_argument("coupled-mode field profile inputs are inconsistent");
    }
    for (std::size_t profile = 0; profile < profiles.size(); ++profile)
    {
        if (profiles[profile]->modes.empty() ||
            profiles[profile]->wavenumbers.size() != profiles[profile]->modes.size() ||
            profiles[profile]->depths.size() < 2 ||
            (profile > 0 && parameters.profileRangesKm[profile] <=
                                  parameters.profileRangesKm[profile - 1]))
        {
            throw std::invalid_argument("coupled-mode MOD profiles are inconsistent");
        }
    }

    PressureField result;
    result.title = parameters.title.empty() ? modes.title : parameters.title;
    result.frequency = modes.frequency;
    result.sourceDepths = parameters.sourceDepths;
    result.receiverDepths = parameters.receiverDepths;
    result.rangesMetres = parameters.rangesMetres;
    result.values.resize(result.sourceDepths.size() * result.receiverDepths.size() *
                         result.rangesMetres.size());
    result.horizontalValues.resize(result.values.size());
    result.verticalValues.resize(result.values.size());
    result.workerThreadsUsed = static_cast<int>(std::min(
        result.sourceDepths.size(),
        static_cast<std::size_t>(std::max(1, parameters.threadCount))));

    constexpr float piSingle = 3.1415926f;
    const ModComplex imaginary(0.0f, 1.0f);
    ModComplex factor = std::sqrt(2.0f * piSingle) *
                        std::exp(imaginary * (piSingle / 4.0f));
    if (parameters.sourceType != 'X')
    {
        factor *= imaginary;
    }
    std::vector<double> interfaces(profiles.size(), 0.0);
    for (std::size_t profile = 1; profile < profiles.size(); ++profile)
    {
        interfaces[profile] = 500.0 *
            (parameters.profileRangesKm[profile - 1] +
             parameters.profileRangesKm[profile]);
    }

    std::vector<std::vector<ModComplex>> receiverModes(profiles.size());
    std::vector<std::vector<ModComplex>> receiverDerivatives(profiles.size());
    for (std::size_t profile = 0; profile < profiles.size(); ++profile)
    {
        const std::size_t modeCount = profiles[profile]->modes.size();
        receiverModes[profile].resize(result.receiverDepths.size() * modeCount);
        receiverDerivatives[profile].resize(result.receiverDepths.size() * modeCount);
        for (std::size_t receiver = 0; receiver < result.receiverDepths.size(); ++receiver)
        {
            const float depth = static_cast<float>(result.receiverDepths[receiver]);
            for (std::size_t mode = 0; mode < modeCount; ++mode)
            {
                const std::size_t index = receiver * modeCount + mode;
                receiverModes[profile][index] =
                    modeAtDepth(*profiles[profile], mode, depth, modes.frequency);
                receiverDerivatives[profile][index] =
                    modeDerivativeAtDepth(*profiles[profile], mode, depth, modes.frequency);
            }
        }
    }

    parallelFor(result.sourceDepths.size(),
                static_cast<std::size_t>(std::max(1, parameters.threadCount)),
                [&](std::size_t sourceIndex) {
        const float sourceSoundSpeed = static_cast<float>(
            sourceIndex < parameters.sourceSoundSpeeds.size()
                ? parameters.sourceSoundSpeeds[sourceIndex] : 1500.0);
        const float omega = 2.0f * piSingle * static_cast<float>(modes.frequency);
        const std::size_t sourceModes = std::min(
            profiles[0]->modes.size(),
            static_cast<std::size_t>(std::max(0, parameters.modeLimit)));
        std::vector<ModComplex> amplitudes(sourceModes);
        for (std::size_t mode = 0; mode < sourceModes; ++mode)
        {
            ModComplex source = modeAtDepth(*profiles[0], mode,
                static_cast<float>(result.sourceDepths[sourceIndex]), modes.frequency);
            if (parameters.beamPattern && sourceIndex == 0)
            {
                const ModComplex k = modValue(profiles[0]->wavenumbers[mode]);
                const float kz2 = std::max(0.0f,
                    std::real(omega * omega / (1500.0f * 1500.0f) - k * k));
                const double angle = (180.0f / piSingle) *
                    std::atan(std::sqrt(kz2) / std::real(k));
                if (!parameters.beamAnglesDegrees.empty())
                {
                    const auto upper = std::upper_bound(
                        parameters.beamAnglesDegrees.begin(),
                        parameters.beamAnglesDegrees.end(), angle);
                    const std::size_t right = std::clamp<std::size_t>(
                        static_cast<std::size_t>(upper - parameters.beamAnglesDegrees.begin()),
                        1, parameters.beamAnglesDegrees.size() - 1);
                    const std::size_t left = right - 1;
                    const double weight = (angle - parameters.beamAnglesDegrees[left]) /
                        (parameters.beamAnglesDegrees[right] -
                         parameters.beamAnglesDegrees[left]);
                    source *= static_cast<float>(
                        (1.0 - weight) * parameters.beamAmplitudes[left] +
                        weight * parameters.beamAmplitudes[right]);
                }
            }
            const ModComplex k = modValue(profiles[0]->wavenumbers[mode]);
            amplitudes[mode] = parameters.sourceType == 'X'
                ? factor * source / k : factor * source / std::sqrt(k);
        }

        std::size_t profile = 0;
        for (std::size_t range = 0; range < result.rangesMetres.size(); ++range)
        {
            const double receiverRange = result.rangesMetres[range];
            const double previousRange = range > 0 ? result.rangesMetres[range - 1] : 0.0;
            if (profile + 1 < profiles.size() &&
                receiverRange > interfaces[profile + 1])
            {
                ++profile;
                advancePhase(amplitudes, *profiles[profile - 1],
                             interfaces[profile] - previousRange);
                amplitudes = projectToNewProfile(
                    *profiles[profile - 1], *profiles[profile], amplitudes,
                    modes.frequency);
                while (profile + 1 < profiles.size() &&
                       receiverRange > interfaces[profile + 1])
                {
                    ++profile;
                    advancePhase(amplitudes, *profiles[profile - 1],
                                 interfaces[profile] - interfaces[profile - 1]);
                    amplitudes = projectToNewProfile(
                        *profiles[profile - 1], *profiles[profile], amplitudes,
                        modes.frequency);
                }
                advancePhase(amplitudes, *profiles[profile],
                             receiverRange - interfaces[profile]);
            }
            else
            {
                advancePhase(amplitudes, *profiles[profile],
                             receiverRange - previousRange);
            }

            for (std::size_t receiver = 0; receiver < result.receiverDepths.size(); ++receiver)
            {
                ModComplex pressure(0.0f, 0.0f);
                ModComplex horizontal(0.0f, 0.0f);
                ModComplex vertical(0.0f, 0.0f);
                const std::size_t profileModeCount = profiles[profile]->modes.size();
                const std::size_t receiverOffset = receiver * profileModeCount;
                for (std::size_t mode = 0; mode < amplitudes.size(); ++mode)
                {
                    const ModComplex receiverMode =
                        receiverModes[profile][receiverOffset + mode];
                    const ModComplex term = amplitudes[mode] * receiverMode;
                    pressure += term;
                    horizontal += sourceSoundSpeed *
                        modValue(profiles[profile]->wavenumbers[mode]) * term / omega;
                    vertical += sourceSoundSpeed * amplitudes[mode] *
                        receiverDerivatives[profile][receiverOffset + mode] /
                        (omega * imaginary);
                }
                if (parameters.sourceType == 'R' && receiverRange != 0.0)
                {
                    pressure /= std::sqrt(static_cast<float>(receiverRange));
                    horizontal /= std::sqrt(static_cast<float>(receiverRange));
                    vertical /= std::sqrt(static_cast<float>(receiverRange));
                }
                const std::size_t index =
                    (sourceIndex * result.receiverDepths.size() + receiver) *
                        result.rangesMetres.size() + range;
                result.values[index] =
                    {pressure.real(), pressure.imag()};
                result.horizontalValues[index] =
                    {horizontal.real(), horizontal.imag()};
                result.verticalValues[index] =
                    {vertical.real(), vertical.imag()};
            }
        }
    });
    return result;
}
}
