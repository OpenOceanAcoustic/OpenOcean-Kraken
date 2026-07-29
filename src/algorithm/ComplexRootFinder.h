#ifndef OPEN_OCEAN_KRAKENC_COMPLEX_ROOT_FINDER_H
#define OPEN_OCEAN_KRAKENC_COMPLEX_ROOT_FINDER_H

#include <complex>
#include <functional>
#include <limits>
#include <vector>

namespace OpenOceanKrakenc
{
struct ScaledComplex
{
    std::complex<double> value{};
    int power10 = 0;
};

enum class RootFailure
{
    None,
    InvalidTolerance,
    NonFiniteValue,
    DegenerateSecant,
    IterationLimit
};

struct RootResult
{
    std::complex<double> root{};
    double log10Residual = 0.0;
    double relativeCorrection = 0.0;
    double absoluteCorrection =
        std::numeric_limits<double>::infinity();
    int iterations = 0;
    bool converged = false;
    RootFailure failure = RootFailure::None;
    std::complex<double> initialGuess{};
};

struct RootConvergenceSpec
{
    double geometricTolerance = 0.0;
    double maximumLog10BackwardError = 0.0;
    double denominatorRelativeTolerance = 0.0;
};

using ScaledFunction = std::function<ScaledComplex(std::complex<double>)>;

class ComplexDispersion;

RootResult complexSecant(
    std::complex<double> initialGuess,
    double tolerance,
    int maxIterations,
    const ScaledFunction &function,
    const RootConvergenceSpec &spec);
RootResult complexSecant(std::complex<double> initialGuess,
                         double tolerance,
                         int maxIterations,
                         const ScaledFunction &function);
RootResult complexSecantDispersion(
    std::complex<double> initialGuess,
    double tolerance,
    int maxIterations,
    const ComplexDispersion &dispersion,
    const std::vector<std::complex<double>> &acceptedRoots,
    const RootConvergenceSpec &spec);
RootResult complexSecantDispersion(
    std::complex<double> initialGuess,
    double tolerance,
    int maxIterations,
    const ComplexDispersion &dispersion,
    const std::vector<std::complex<double>> &acceptedRoots);
}

#endif
