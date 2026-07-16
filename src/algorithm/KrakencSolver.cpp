#include "algorithm/KrakencSolver.h"

#include "algorithm/ComplexDispersion.h"
#include "algorithm/ComplexMatrixBuilder.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace OpenOceanKrakenc
{
namespace
{
double halton(std::size_t index, std::size_t base)
{
    double fraction = 1.0;
    double result = 0.0;
    while (index > 0)
    {
        fraction /= static_cast<double>(base);
        result += fraction * static_cast<double>(index % base);
        index /= base;
    }
    return result;
}

bool finite(std::complex<double> value)
{
    return std::isfinite(value.real()) && std::isfinite(value.imag());
}

double minimumSoundSpeed(const AcousticCase &input)
{
    double result = std::numeric_limits<double>::infinity();
    for (const AcousticLayer &layer : input.layers)
    {
        for (const AcousticSample &sample : layer.samples)
        {
            result = std::min(result, sample.cp);
        }
    }
    if (input.bottom.type == AcousticBoundaryType::HalfSpace)
    {
        result = std::min(result, input.bottom.cp);
    }
    if (input.top.type == AcousticBoundaryType::HalfSpace)
    {
        result = std::min(result, input.top.cp);
    }
    return result;
}

bool duplicateRoot(const std::vector<std::complex<double>> &roots,
                   std::complex<double> candidate,
                   double tolerance)
{
    return std::any_of(roots.begin(), roots.end(), [&](std::complex<double> root) {
        return std::abs(root - candidate) <=
               tolerance * std::max({1.0, std::abs(root), std::abs(candidate)});
    });
}

std::vector<ModeRoot> solveMesh(const AcousticCase &input,
                                const AcousticMatrix &matrix,
                                const std::vector<std::complex<double>> &initialGuesses = {})
{
    const ComplexDispersion dispersion(input, matrix);
    const double adjustedCLow = std::max(input.cLow, 0.99 * minimumSoundSpeed(input));
    const double upperK = matrix.omega / adjustedCLow;
    const double lowerK = matrix.omega / input.cHigh;
    const double totalDepth = input.layers.back().bottomDepth -
                              input.layers.front().topDepth;
    const std::size_t maxModes = initialGuesses.empty()
                                     ? std::max<std::size_t>(
                                           16, static_cast<std::size_t>(
                                                   std::ceil(2.0 * totalDepth * input.frequency /
                                                             adjustedCLow * 1.1)))
                                     : initialGuesses.size();
    const int maxRestarts = input.enableRootRestarts
                                ? std::max(10, static_cast<int>(matrix.b1.size() / 50))
                                : 0;

    std::vector<std::complex<double>> accepted;
    std::vector<ModeRoot> modes;
    std::complex<double> guess(upperK * upperK, 0.0);
    int restartCount = 0;

    while (modes.size() < maxModes)
    {
        if (modes.size() < initialGuesses.size())
        {
            guess = initialGuesses[modes.size()];
        }
        else
        {
            guess *= 1.00001;
        }
        const double tolerance = std::max(
            1.0e-14, std::abs(guess) * static_cast<double>(matrix.b1.size()) * 1.0e-14);
        RootResult diagnostic = complexSecantDispersion(
            guess, tolerance, 1000, dispersion, accepted);
        const std::complex<double> candidate = diagnostic.root;
        const std::complex<double> k = finite(candidate)
                                           ? std::sqrt(candidate)
                                           : std::complex<double>(
                                                 std::numeric_limits<double>::quiet_NaN(), 0.0);
        const bool valid = diagnostic.converged && finite(candidate) && finite(k) &&
                           k.real() >= lowerK && k.real() <= upperK * 1.01 &&
                           !duplicateRoot(accepted, candidate, 1.0e-8);

        if (!valid)
        {
            if (restartCount >= maxRestarts)
            {
                break;
            }
            ++restartCount;
            const auto point = deterministicRestartPoint(static_cast<std::size_t>(restartCount));
            const std::complex<double> kz(
                point.first * upperK, 0.01 * point.second * upperK);
            guess = upperK * upperK - kz * kz;
            continue;
        }

        ModeRoot mode;
        mode.eigenvalue = candidate;
        mode.wavenumber = k;
        mode.diagnostic = diagnostic;
        mode.discoveryOrder = modes.size();
        accepted.push_back(candidate);
        modes.push_back(mode);
        guess = candidate;
    }

    std::stable_sort(modes.begin(), modes.end(), [](const ModeRoot &left, const ModeRoot &right) {
        if (left.eigenvalue.real() != right.eigenvalue.real())
        {
            return left.eigenvalue.real() > right.eigenvalue.real();
        }
        if (left.eigenvalue.imag() != right.eigenvalue.imag())
        {
            return left.eigenvalue.imag() > right.eigenvalue.imag();
        }
        return left.discoveryOrder < right.discoveryOrder;
    });
    return modes;
}

std::complex<double> extrapolatedInitialGuess(
    const std::vector<int> &multipliers,
    const std::vector<std::complex<double>> &values,
    int targetMultiplier)
{
    std::vector<std::complex<double>> work = values;
    const double target = 1.0 /
                          static_cast<double>(targetMultiplier * targetMultiplier);
    for (std::size_t order = 1; order < work.size(); ++order)
    {
        for (std::size_t index = 0; index + order < work.size(); ++index)
        {
            const double x1 = 1.0 / static_cast<double>(
                                        multipliers[index] * multipliers[index]);
            const double x2 = 1.0 / static_cast<double>(
                                        multipliers[index + order] *
                                        multipliers[index + order]);
            work[index] = ((target - x2) * work[index] -
                           (target - x1) * work[index + 1]) /
                          (x1 - x2);
        }
    }
    return work.front();
}
}

std::pair<double, double> deterministicRestartPoint(std::size_t index)
{
    if (index == 0)
    {
        throw std::invalid_argument("deterministic restart index is one-based");
    }
    return {halton(index, 2), halton(index, 3)};
}

std::complex<double> richardsonExtrapolate(
    const std::vector<int> &multipliers,
    const std::vector<std::complex<double>> &values)
{
    if (multipliers.empty() || multipliers.size() != values.size())
    {
        throw std::invalid_argument("Richardson inputs must be nonempty and equally sized");
    }
    std::vector<double> x;
    x.reserve(multipliers.size());
    for (int multiplier : multipliers)
    {
        if (multiplier < 1)
        {
            throw std::invalid_argument("Richardson multipliers must be positive");
        }
        x.push_back(1.0 / static_cast<double>(multiplier * multiplier));
    }

    std::vector<std::complex<double>> work = values;
    for (std::size_t order = 1; order < work.size(); ++order)
    {
        for (std::size_t i = 0; i + order < work.size(); ++i)
        {
            const double denominator = x[i] - x[i + order];
            if (denominator == 0.0)
            {
                throw std::invalid_argument("Richardson multipliers must be unique");
            }
            work[i] = (-x[i + order] * work[i] + x[i] * work[i + 1]) /
                      denominator;
        }
    }
    return work.front();
}

AcousticSolveResult solveAcousticModes(const AcousticCase &input)
{
    static const std::vector<int> meshMultipliers{1, 2, 4, 8, 16};
    AcousticSolveResult result;
    std::vector<std::vector<ModeRoot>> meshModes;

    for (int multiplier : meshMultipliers)
    {
        const AcousticMatrix matrix = buildAcousticMatrix(input, multiplier);
        std::vector<std::complex<double>> initialGuesses;
        if (meshModes.size() >= 2)
        {
            const std::size_t modeCount = meshModes.back().size();
            initialGuesses.reserve(modeCount);
            const std::vector<int> priorMultipliers(
                meshMultipliers.begin(),
                meshMultipliers.begin() + static_cast<std::ptrdiff_t>(meshModes.size()));
            for (std::size_t mode = 0; mode < modeCount; ++mode)
            {
                std::vector<std::complex<double>> values;
                values.reserve(meshModes.size());
                bool complete = true;
                for (const auto &set : meshModes)
                {
                    if (mode >= set.size())
                    {
                        complete = false;
                        break;
                    }
                    values.push_back(set[mode].eigenvalue);
                }
                if (!complete)
                {
                    break;
                }
                initialGuesses.push_back(extrapolatedInitialGuess(
                    priorMultipliers, values, multiplier));
            }
        }
        std::vector<ModeRoot> current = solveMesh(input, matrix, initialGuesses);
        if (current.empty())
        {
            throw std::runtime_error("complex acoustic root search found no modes");
        }
        meshModes.push_back(std::move(current));
        if (meshModes.size() == 2 &&
            meshModes.front().size() < meshModes.back().size())
        {
            std::vector<std::complex<double>> reverseGuesses;
            reverseGuesses.reserve(meshModes.back().size());
            for (const ModeRoot &mode : meshModes.back())
            {
                reverseGuesses.push_back(mode.eigenvalue);
            }
            const AcousticMatrix coarseMatrix = buildAcousticMatrix(
                input, meshMultipliers.front());
            std::vector<ModeRoot> recovered = solveMesh(
                input, coarseMatrix, reverseGuesses);
            if (recovered.size() > meshModes.front().size())
            {
                meshModes.front() = std::move(recovered);
                result.modes = meshModes.front();
                result.meshModeCounts.front() =
                    static_cast<int>(meshModes.front().size());
                result.meshWavenumberSets.front().clear();
                for (const ModeRoot &mode : meshModes.front())
                {
                    result.meshWavenumberSets.front().push_back(mode.wavenumber);
                }
            }
        }
        result.meshModeCounts.push_back(
            static_cast<int>(meshModes.back().size()));
        result.meshWavenumberSets.emplace_back();
        for (const ModeRoot &mode : meshModes.back())
        {
            result.meshWavenumberSets.back().push_back(mode.wavenumber);
        }
        result.meshSetsUsed = static_cast<int>(meshModes.size());

        if (meshModes.size() == 1)
        {
            result.modes = meshModes.front();
            result.meshErrors.push_back(std::numeric_limits<double>::infinity());
            if (input.rMaxKm == 0.0)
            {
                break;
            }
            continue;
        }

        std::size_t commonModes = meshModes.front().size();
        for (const auto &set : meshModes)
        {
            commonModes = std::min(commonModes, set.size());
        }
        if (commonModes == 0)
        {
            throw std::runtime_error("mesh sets have no common acoustic modes");
        }

        std::vector<ModeRoot> extrapolated;
        extrapolated.reserve(commonModes);
        const std::vector<int> usedMultipliers(
            meshMultipliers.begin(), meshMultipliers.begin() + result.meshSetsUsed);
        for (std::size_t mode = 0; mode < commonModes; ++mode)
        {
            std::vector<std::complex<double>> values;
            values.reserve(meshModes.size());
            for (const auto &set : meshModes)
            {
                values.push_back(set[mode].eigenvalue);
            }
            ModeRoot output = meshModes.back()[mode];
            output.eigenvalue = richardsonExtrapolate(usedMultipliers, values);
            output.wavenumber = std::sqrt(output.eigenvalue);
            extrapolated.push_back(output);
        }

        const std::size_t key = std::min(commonModes - 1, 2 * commonModes / 3);
        const double error = std::abs(extrapolated[key].eigenvalue -
                                      result.modes[std::min(key, result.modes.size() - 1)].eigenvalue);
        result.meshErrors.push_back(error);
        result.modes = std::move(extrapolated);
        if (error * 1000.0 * input.rMaxKm < 1.0)
        {
            break;
        }
    }
    return result;
}
}
