#ifndef OPEN_OCEAN_KRAKENC_SOLVER_H
#define OPEN_OCEAN_KRAKENC_SOLVER_H

#include "algorithm/AcousticCase.h"
#include "algorithm/ComplexRootFinder.h"

#include <complex>
#include <cstddef>
#include <utility>
#include <vector>

namespace OpenOceanKrakenc
{
struct ModeRoot
{
    std::complex<double> eigenvalue{};
    std::complex<double> wavenumber{};
    RootResult diagnostic;
    std::size_t discoveryOrder = 0;
};

struct AcousticSolveResult
{
    std::vector<ModeRoot> modes;
    std::vector<double> meshErrors;
    std::vector<int> meshModeCounts;
    std::vector<std::vector<std::complex<double>>> meshWavenumberSets;
    int meshSetsUsed = 0;
};

std::pair<double, double> deterministicRestartPoint(std::size_t index);

std::complex<double> richardsonExtrapolate(
    const std::vector<int> &multipliers,
    const std::vector<std::complex<double>> &values);

AcousticSolveResult solveAcousticModes(const AcousticCase &input);
}

#endif
