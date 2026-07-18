#include "algorithm/KrakencSolver.h"

#include "algorithm/ComplexDispersion.h"
#include "algorithm/ComplexMatrixBuilder.h"
#include "algorithm/ComplexModeNormalization.h"
#include "algorithm/ComplexModeSolver.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>
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

bool duplicateRoot(const std::vector<std::complex<double>> &roots,
                   std::complex<double> candidate,
                   double tolerance)
{
    return std::any_of(roots.begin(), roots.end(), [&](std::complex<double> root) {
        return std::abs(root - candidate) <=
               tolerance * std::max({1.0, std::abs(root), std::abs(candidate)});
    });
}

struct MeshSolveResult
{
    std::vector<ModeRoot> modes;
    std::vector<RootSearchDiagnostic> diagnostics;
};

MeshSolveResult solveMesh(
    const AcousticCase &input,
    const AcousticMatrix &matrix,
    int meshMultiplier,
    const std::vector<std::complex<double>> &initialGuesses = {})
{
    const ComplexDispersion dispersion(input, matrix);
    const SpectralMinimumPhaseSpeed spectralFloor =
        spectralMinimumPhaseSpeed(input);
    const double adjustedCLow = std::max(input.cLow, 0.99 * spectralFloor.speed);
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
    const int maxTries = input.enableRootRestarts && meshMultiplier <= 2
                             ? std::max(10, static_cast<int>(matrix.b1.size() / 50))
                             : 1;
    const int maxRestarts = maxTries - 1;

    std::vector<std::complex<double>> accepted;
    MeshSolveResult result;
    std::complex<double> guess(upperK * upperK, 0.0);
    int restartCount = 0;

    while (result.modes.size() < maxModes)
    {
        if (result.modes.size() < initialGuesses.size())
        {
            guess = initialGuesses[result.modes.size()];
        }
        else
        {
            const bool denseElasticStack =
                spectralFloor.elasticPresent && input.layers.size() >= 4;
            // The frozen Fortran build rounds the unkinded step literal to
            // REAL before promotion. That rounding selects the compatible
            // branch in dense elastic stacks; simpler stacks retain their
            // established C++ root sequence and regression baseline.
            guess *= denseElasticStack
                         ? static_cast<double>(1.00001F)
                         : 1.00001;
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

        RootSearchDiagnostic audit;
        audit.meshMultiplier = meshMultiplier;
        audit.initialGuess = diagnostic.initialGuess;
        audit.candidate = candidate;
        audit.log10Residual = diagnostic.log10Residual;
        audit.relativeCorrection = diagnostic.relativeCorrection;
        audit.searchLowerWavenumber = lowerK;
        audit.searchUpperWavenumber = upperK;
        audit.accepted = valid;
        if (!diagnostic.converged)
        {
            audit.decision = "rejected: secant did not converge";
        }
        else if (!finite(candidate) || !finite(k))
        {
            audit.decision = "rejected: non-finite root";
        }
        else if (k.real() < lowerK || k.real() > upperK * 1.01)
        {
            audit.decision = "rejected: outside search interval";
        }
        else if (duplicateRoot(accepted, candidate, 1.0e-8))
        {
            audit.decision = "rejected: duplicate root";
        }
        else
        {
            audit.decision = "accepted: converged, finite, in-band and unique";
        }
        result.diagnostics.push_back(audit);

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
        mode.discoveryOrder = result.modes.size();
        mode.meshMultiplier = meshMultiplier;
        mode.searchLowerWavenumber = lowerK;
        mode.searchUpperWavenumber = upperK;
        mode.acceptanceDiagnostic = audit.decision;
        accepted.push_back(candidate);
        result.modes.push_back(mode);
        guess = candidate;
    }

    std::stable_sort(result.modes.begin(), result.modes.end(), [](const ModeRoot &left, const ModeRoot &right) {
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
    return result;
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

using ModeTrack = std::vector<const ModeRoot *>;

double relativeWavenumberDistance(const ModeRoot &left,
                                  const ModeRoot &right)
{
    return std::abs(left.wavenumber - right.wavenumber) /
           std::max({std::abs(left.wavenumber),
                     std::abs(right.wavenumber),
                     std::numeric_limits<double>::min()});
}

// Retained as an audit helper for physical nearest-root continuity; the
// KrakenC-compatible solve path below intentionally uses ordinal tracks.
[[maybe_unused]] std::vector<ModeTrack> matchedModeTracks(
    const std::vector<std::vector<ModeRoot>> &meshModes)
{
    if (meshModes.empty())
    {
        return {};
    }
    std::vector<ModeTrack> tracks;
    tracks.reserve(meshModes.front().size());
    for (const ModeRoot &mode : meshModes.front())
    {
        tracks.push_back({&mode});
    }

    constexpr double gapPenalty = 0.01;
    constexpr double maximumMatchDistance = 0.05;
    for (std::size_t setIndex = 1; setIndex < meshModes.size(); ++setIndex)
    {
        const std::vector<ModeRoot> &current = meshModes[setIndex];
        const std::size_t rows = tracks.size() + 1;
        const std::size_t columns = current.size() + 1;
        const double infinity = std::numeric_limits<double>::infinity();
        std::vector<double> cost(rows * columns, infinity);
        std::vector<char> action(rows * columns, 0);
        const auto index = [columns](std::size_t row, std::size_t column) {
            return row * columns + column;
        };
        cost[index(0, 0)] = 0.0;
        for (std::size_t row = 1; row < rows; ++row)
        {
            cost[index(row, 0)] = cost[index(row - 1, 0)] + gapPenalty;
            action[index(row, 0)] = 'T';
        }
        for (std::size_t column = 1; column < columns; ++column)
        {
            cost[index(0, column)] = cost[index(0, column - 1)] + gapPenalty;
            action[index(0, column)] = 'C';
        }
        for (std::size_t row = 1; row < rows; ++row)
        {
            for (std::size_t column = 1; column < columns; ++column)
            {
                double best = cost[index(row - 1, column)] + gapPenalty;
                char bestAction = 'T';
                const double skipCurrent =
                    cost[index(row, column - 1)] + gapPenalty;
                if (skipCurrent < best)
                {
                    best = skipCurrent;
                    bestAction = 'C';
                }
                const double distance = relativeWavenumberDistance(
                    *tracks[row - 1].back(), current[column - 1]);
                const double match = cost[index(row - 1, column - 1)] + distance;
                if (distance <= maximumMatchDistance && match <= best)
                {
                    best = match;
                    bestAction = 'M';
                }
                cost[index(row, column)] = best;
                action[index(row, column)] = bestAction;
            }
        }

        std::vector<std::pair<std::size_t, std::size_t>> pairs;
        std::size_t row = tracks.size();
        std::size_t column = current.size();
        while (row > 0 || column > 0)
        {
            const char step = action[index(row, column)];
            if (step == 'M')
            {
                pairs.emplace_back(row - 1, column - 1);
                --row;
                --column;
            }
            else if (step == 'T')
            {
                --row;
            }
            else
            {
                --column;
            }
        }
        std::reverse(pairs.begin(), pairs.end());
        std::vector<ModeTrack> next;
        next.reserve(pairs.size());
        for (const auto [trackIndex, modeIndex] : pairs)
        {
            ModeTrack track = tracks[trackIndex];
            track.push_back(&current[modeIndex]);
            next.push_back(std::move(track));
        }
        tracks = std::move(next);
    }
    return tracks;
}

std::vector<ModeTrack> ordinalModeTracks(
    const std::vector<std::vector<ModeRoot>> &meshModes)
{
    if (meshModes.empty())
    {
        return {};
    }
    std::size_t commonModes = meshModes.front().size();
    for (const std::vector<ModeRoot> &modes : meshModes)
    {
        commonModes = std::min(commonModes, modes.size());
    }
    std::vector<ModeTrack> tracks(commonModes);
    for (std::size_t modeIndex = 0; modeIndex < commonModes; ++modeIndex)
    {
        tracks[modeIndex].reserve(meshModes.size());
        for (const std::vector<ModeRoot> &modes : meshModes)
        {
            tracks[modeIndex].push_back(&modes[modeIndex]);
        }
    }
    return tracks;
}

double trackRelativeSpread(const ModeTrack &track)
{
    double result = 0.0;
    for (std::size_t left = 0; left < track.size(); ++left)
    {
        for (std::size_t right = left + 1; right < track.size(); ++right)
        {
            result = std::max(
                result,
                relativeWavenumberDistance(*track[left], *track[right]));
        }
    }
    return result;
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
            const std::vector<ModeTrack> tracks = ordinalModeTracks(meshModes);
            initialGuesses.reserve(tracks.size());
            const std::vector<int> priorMultipliers(
                meshMultipliers.begin(),
                meshMultipliers.begin() + static_cast<std::ptrdiff_t>(meshModes.size()));
            for (const ModeTrack &track : tracks)
            {
                std::vector<std::complex<double>> values;
                values.reserve(track.size());
                for (const ModeRoot *mode : track)
                {
                    values.push_back(mode->eigenvalue);
                }
                initialGuesses.push_back(extrapolatedInitialGuess(
                    priorMultipliers, values, multiplier));
            }
        }
        MeshSolveResult meshSolve = solveMesh(
            input, matrix, multiplier, initialGuesses);
        if (meshSolve.modes.empty())
        {
            throw std::runtime_error("complex acoustic root search found no modes");
        }
        result.rootSearchDiagnostics.insert(
            result.rootSearchDiagnostics.end(),
            meshSolve.diagnostics.begin(), meshSolve.diagnostics.end());
        meshModes.push_back(std::move(meshSolve.modes));
        if (meshModes.size() == 1)
        {
            for (std::size_t modeIndex = 0;
                 modeIndex < meshModes.front().size(); ++modeIndex)
            {
                ModeRoot &mode = meshModes.front()[modeIndex];
                mode.coarseModeIndex = modeIndex;
                mode.coarseEigenvalue = mode.eigenvalue;
                mode.hasCoarseProvenance = true;
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

        const std::vector<ModeTrack> tracks = ordinalModeTracks(meshModes);
        if (tracks.empty())
        {
            throw std::runtime_error("mesh sets have no common acoustic modes");
        }

        std::vector<ModeRoot> extrapolated;
        extrapolated.reserve(tracks.size());
        const std::vector<int> usedMultipliers(
            meshMultipliers.begin(), meshMultipliers.begin() + result.meshSetsUsed);
        for (const ModeTrack &track : tracks)
        {
            std::vector<std::complex<double>> values;
            values.reserve(track.size());
            for (const ModeRoot *mode : track)
            {
                values.push_back(mode->eigenvalue);
            }
            ModeRoot output = *track.back();
            output.eigenvalue = richardsonExtrapolate(usedMultipliers, values);
            output.wavenumber = std::sqrt(output.eigenvalue);
            output.meshRelativeSpread = trackRelativeSpread(track);
            output.matchedMeshSets = static_cast<int>(track.size());
            output.coarseModeIndex = static_cast<std::size_t>(
                track.front() - meshModes.front().data());
            output.coarseEigenvalue = track.front()->eigenvalue;
            output.hasCoarseProvenance = true;
            std::ostringstream diagnostic;
            diagnostic << "accepted: Fortran ordinal correspondence across "
                       << track.size() << " mesh sets; relative spread "
                       << output.meshRelativeSpread;
            output.acceptanceDiagnostic = diagnostic.str();
            extrapolated.push_back(output);
        }

        const std::size_t key = std::min(
            extrapolated.size() - 1, 2 * extrapolated.size() / 3);
        const double error = std::abs(extrapolated[key].eigenvalue -
                                      result.modes[std::min(key, result.modes.size() - 1)].eigenvalue);
        result.meshErrors.push_back(error);
        result.modes = std::move(extrapolated);
        if (error * 1000.0 * input.rMaxKm < 1.0)
        {
            break;
        }
    }
    const bool hasRoughness = input.bottom.roughnessRms != 0.0 ||
        std::any_of(input.layers.begin(), input.layers.end(),
                    [](const AcousticLayer &layer) {
                        return layer.roughnessRms != 0.0;
                    });
    if (hasRoughness)
    {
        const AcousticMatrix matrix = buildAcousticMatrix(input, 1);
        for (ModeRoot &root : result.modes)
        {
            const ComplexModeResult raw = solveAcousticMode(
                input, matrix, root.eigenvalue);
            if (!raw.converged)
            {
                throw std::runtime_error(
                    "roughness mode extraction failed before scatter correction");
            }
            const NormalizedComplexMode normalized = normalizeAcousticMode(
                input, matrix, root.eigenvalue, raw.turningPoint, raw.mode);
            root.scatterPerturbation = normalized.scatterPerturbation;
            root.wavenumber = std::sqrt(
                root.eigenvalue + root.scatterPerturbation);
        }
    }
    for (ModeRoot &root : result.modes)
    {
        if (root.wavenumber.imag() > 0.0)
        {
            root.wavenumber = {root.wavenumber.real(), 0.0};
        }
    }
    return result;
}

SpectralMinimumPhaseSpeed spectralMinimumPhaseSpeed(
    const AcousticCase &input)
{
    SpectralMinimumPhaseSpeed result;
    result.speed = std::numeric_limits<double>::infinity();
    for (const AcousticLayer &layer : input.layers)
    {
        for (const AcousticSample &sample : layer.samples)
        {
            if (sample.cs > 0.0)
            {
                result.elasticPresent = true;
                result.speed = std::min(result.speed, sample.cs);
            }
            else
            {
                result.speed = std::min(result.speed, sample.cp);
            }
        }
    }
    const auto includeHalfSpace = [&](const AcousticBoundary &boundary) {
        if (boundary.type != AcousticBoundaryType::HalfSpace)
        {
            return;
        }
        if (boundary.cs > 0.0)
        {
            result.elasticPresent = true;
            result.speed = std::min(result.speed, boundary.cs);
        }
        else
        {
            result.speed = std::min(result.speed, boundary.cp);
        }
    };
    includeHalfSpace(input.bottom);
    includeHalfSpace(input.top);
    if (result.elasticPresent)
    {
        result.speed *= 0.85;
    }
    if (!std::isfinite(result.speed) || !(result.speed > 0.0))
    {
        throw std::invalid_argument(
            "cannot derive a positive spectral phase-speed floor");
    }
    return result;
}
}
