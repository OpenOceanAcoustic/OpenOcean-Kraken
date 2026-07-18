#ifndef OPEN_OCEAN_KRAKENC_SOLVER_H
#define OPEN_OCEAN_KRAKENC_SOLVER_H

#include "algorithm/AcousticCase.h"
#include "algorithm/ComplexRootFinder.h"

#include <complex>
#include <cstddef>
#include <string>
#include <utility>
#include <vector>

namespace OpenOceanKrakenc
{
struct SpectralMinimumPhaseSpeed
{
    double speed = 0.0;
    bool elasticPresent = false;
};

struct ModeRoot
{
    std::complex<double> eigenvalue{};
    std::complex<double> wavenumber{};
    std::complex<double> scatterPerturbation{};
    RootResult diagnostic;
    std::size_t discoveryOrder = 0;
    int meshMultiplier = 1;
    double searchLowerWavenumber = 0.0;
    double searchUpperWavenumber = 0.0;
    double meshRelativeSpread = 0.0;
    int matchedMeshSets = 1;
    std::size_t coarseModeIndex = 0;
    std::complex<double> coarseEigenvalue{};
    bool hasCoarseProvenance = false;
    std::string acceptanceDiagnostic;
};

struct RootSearchDiagnostic
{
    int meshMultiplier = 1;
    std::complex<double> initialGuess{};
    std::complex<double> candidate{};
    double log10Residual = 0.0;
    double relativeCorrection = 0.0;
    double searchLowerWavenumber = 0.0;
    double searchUpperWavenumber = 0.0;
    bool accepted = false;
    std::string decision;
};

struct AcousticSolveResult
{
    std::vector<ModeRoot> modes;
    std::vector<double> meshErrors;
    std::vector<int> meshModeCounts;
    std::vector<std::vector<std::complex<double>>> meshWavenumberSets;
    std::vector<RootSearchDiagnostic> rootSearchDiagnostics;
    int meshSetsUsed = 0;
};

std::pair<double, double> deterministicRestartPoint(std::size_t index);

std::complex<double> richardsonExtrapolate(
    const std::vector<int> &multipliers,
    const std::vector<std::complex<double>> &values);

SpectralMinimumPhaseSpeed spectralMinimumPhaseSpeed(
    const AcousticCase &input);

AcousticSolveResult solveAcousticModes(const AcousticCase &input);
}

#endif
