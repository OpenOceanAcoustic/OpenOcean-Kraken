#include "algorithm/ComplexRootFinder.h"

#include "algorithm/ComplexDispersion.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>

namespace OpenOceanKrakenc
{
namespace
{
bool finite(std::complex<double> value)
{
    return std::isfinite(value.real()) && std::isfinite(value.imag());
}

double residualLog10(const ScaledComplex &value)
{
    const double magnitude = std::abs(value.value);
    if (magnitude == 0.0)
    {
        return -std::numeric_limits<double>::infinity();
    }
    if (!std::isfinite(magnitude))
    {
        return std::numeric_limits<double>::infinity();
    }
    return std::log10(magnitude) + static_cast<double>(value.power10);
}

RootResult failureResult(std::complex<double> initialGuess,
                         std::complex<double> root,
                         int iterations,
                         RootFailure failure)
{
    RootResult result;
    result.initialGuess = initialGuess;
    result.root = root;
    result.iterations = iterations;
    result.failure = failure;
    result.log10Residual = std::numeric_limits<double>::infinity();
    result.relativeCorrection = std::numeric_limits<double>::infinity();
    return result;
}

double power10(int exponent)
{
    constexpr int minimum = -308;
    constexpr int maximum = 308;
    static const std::array<double, maximum - minimum + 1> values = [] {
        std::array<double, maximum - minimum + 1> result{};
        for (int value = minimum; value <= maximum; ++value)
        {
            result[static_cast<std::size_t>(value - minimum)] =
                std::pow(10.0, value);
        }
        return result;
    }();
    if (exponent >= minimum && exponent <= maximum)
    {
        return values[static_cast<std::size_t>(exponent - minimum)];
    }
    return std::pow(10.0, exponent);
}
}

template <typename Function>
RootResult complexSecantImpl(std::complex<double> initialGuess,
                             double tolerance,
                             int maxIterations,
                             const Function &function)
{
    if (!(tolerance > 0.0) || !std::isfinite(tolerance) || maxIterations < 1)
    {
        return failureResult(initialGuess, initialGuess, 0,
                             tolerance > 0.0 ? RootFailure::IterationLimit
                                             : RootFailure::InvalidTolerance);
    }

    std::complex<double> x2 = initialGuess;
    std::complex<double> x1 = x2 + 100.0 * tolerance;
    ScaledComplex f1 = function(x1);
    if (!finite(f1.value))
    {
        return failureResult(initialGuess, x2, 0, RootFailure::NonFiniteValue);
    }

    for (int iteration = 1; iteration <= maxIterations; ++iteration)
    {
        const std::complex<double> x0 = x1;
        const ScaledComplex f0 = f1;
        x1 = x2;
        f1 = function(x1);
        if (!finite(f1.value))
        {
            return failureResult(initialGuess, x1, iteration, RootFailure::NonFiniteValue);
        }

        const int exponentDifference = f0.power10 - f1.power10;
        const double scale = power10(exponentDifference);
        const std::complex<double> numerator = f1.value * (x1 - x0);
        const std::complex<double> denominator = f1.value - f0.value * scale;

        std::complex<double> shift;
        if (!finite(numerator) || !finite(denominator))
        {
            return failureResult(initialGuess, x1, iteration, RootFailure::NonFiniteValue);
        }
        if (denominator == std::complex<double>(0.0, 0.0) ||
            std::norm(numerator) >= std::norm(denominator) * std::norm(x1))
        {
            shift = 0.1 * tolerance;
        }
        else
        {
            shift = numerator / denominator;
        }

        x2 = x1 - shift;
        if (!finite(x2))
        {
            return failureResult(initialGuess, x1, iteration, RootFailure::NonFiniteValue);
        }

        if (std::abs(x2 - x1) + std::abs(x2 - x0) < tolerance)
        {
            const double relativeCorrection =
                (std::abs(x2 - x1) + std::abs(x2 - x0)) /
                std::max(std::abs(x2), std::numeric_limits<double>::min());
            const ScaledComplex finalValue = function(x2);
            if (!finite(finalValue.value))
            {
                return failureResult(initialGuess, x2, iteration, RootFailure::NonFiniteValue);
            }
            RootResult result;
            result.initialGuess = initialGuess;
            result.root = x2;
            result.log10Residual = residualLog10(finalValue);
            result.relativeCorrection = relativeCorrection;
            result.iterations = iteration;
            result.converged = true;
            result.failure = RootFailure::None;
            return result;
        }
    }

    RootResult result = failureResult(initialGuess, x2, maxIterations,
                                      RootFailure::IterationLimit);
    const ScaledComplex finalValue = function(x2);
    if (finite(finalValue.value))
    {
        result.log10Residual = residualLog10(finalValue);
    }
    return result;
}

RootResult complexSecant(std::complex<double> initialGuess,
                         double tolerance,
                         int maxIterations,
                         const ScaledFunction &function)
{
    return complexSecantImpl(initialGuess, tolerance, maxIterations, function);
}

RootResult complexSecantDispersion(
    std::complex<double> initialGuess,
    double tolerance,
    int maxIterations,
    const ComplexDispersion &dispersion,
    const std::vector<std::complex<double>> &acceptedRoots)
{
    const auto function = [&](std::complex<double> value) {
        return dispersion.evaluate(value, acceptedRoots);
    };
    return complexSecantImpl(initialGuess, tolerance, maxIterations, function);
}
}
