#include "algorithm/ComplexRootFinder.h"

#include "algorithm/ComplexDispersion.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
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

bool exactZero(const ScaledComplex &value)
{
    return std::abs(value.value) == 0.0;
}

double backwardErrorLog10(const ScaledComplex &value,
                          const ScaledComplex &slope,
                          bool hasSlope)
{
    if (exactZero(value))
    {
        return -std::numeric_limits<double>::infinity();
    }
    const double valueMagnitude = std::abs(value.value);
    const double slopeMagnitude = std::abs(slope.value);
    if (!hasSlope || !(valueMagnitude > 0.0) ||
        !(slopeMagnitude > 0.0) || !std::isfinite(valueMagnitude) ||
        !std::isfinite(slopeMagnitude))
    {
        return std::numeric_limits<double>::infinity();
    }
    return std::log10(valueMagnitude) - std::log10(slopeMagnitude) +
           static_cast<double>(value.power10) -
           static_cast<double>(slope.power10);
}

double magnitudeFromLog10(double logarithm)
{
    if (logarithm == -std::numeric_limits<double>::infinity())
    {
        return 0.0;
    }
    if (!std::isfinite(logarithm) ||
        logarithm > std::log10(std::numeric_limits<double>::max()))
    {
        return std::numeric_limits<double>::infinity();
    }
    return std::pow(10.0, logarithm);
}

bool makeScaledSlope(const std::complex<double> &denominator,
                     int referencePower,
                     const std::complex<double> &delta,
                     ScaledComplex &slope)
{
    const double denominatorMagnitude = std::abs(denominator);
    const double deltaMagnitude = std::abs(delta);
    if (!(denominatorMagnitude > 0.0) || !(deltaMagnitude > 0.0) ||
        !std::isfinite(denominatorMagnitude) ||
        !std::isfinite(deltaMagnitude))
    {
        return false;
    }

    const double log10Mantissa =
        std::log10(denominatorMagnitude) - std::log10(deltaMagnitude);
    const double exponentAsDouble = std::floor(log10Mantissa);
    if (exponentAsDouble <
            static_cast<double>(std::numeric_limits<int>::min()) ||
        exponentAsDouble >
            static_cast<double>(std::numeric_limits<int>::max()))
    {
        return false;
    }
    const int exponent = static_cast<int>(exponentAsDouble);
    const long long combinedPower =
        static_cast<long long>(referencePower) + exponent;
    if (combinedPower < std::numeric_limits<int>::min() ||
        combinedPower > std::numeric_limits<int>::max())
    {
        return false;
    }

    const std::complex<double> direction =
        (denominator / denominatorMagnitude) /
        (delta / deltaMagnitude);
    const double mantissaMagnitude =
        std::pow(10.0, log10Mantissa - exponentAsDouble);
    slope.value = direction * mantissaMagnitude;
    slope.power10 = static_cast<int>(combinedPower);
    return finite(slope.value) && std::abs(slope.value) > 0.0;
}

RootResult makeFailure(std::complex<double> initialGuess,
                       std::complex<double> root,
                       int iterations,
                       RootFailure failure,
                       double log10Residual,
                       double absoluteCorrection)
{
    RootResult result;
    result.initialGuess = initialGuess;
    result.root = root;
    result.iterations = iterations;
    result.failure = failure;
    result.log10Residual = log10Residual;
    result.absoluteCorrection = absoluteCorrection;
    result.relativeCorrection =
        absoluteCorrection /
        std::max(std::abs(root), std::numeric_limits<double>::min());
    return result;
}

RootResult makeSuccess(std::complex<double> initialGuess,
                       std::complex<double> root,
                       int iterations,
                       double log10Residual,
                       double absoluteCorrection)
{
    RootResult result = makeFailure(
        initialGuess, root, iterations, RootFailure::None,
        log10Residual, absoluteCorrection);
    result.converged = true;
    return result;
}

RootConvergenceSpec defaultConvergenceSpec(double tolerance)
{
    return {
        tolerance,
        std::log10(tolerance),
        64.0 * std::numeric_limits<double>::epsilon()};
}

double power10(std::int64_t exponent)
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
    static const std::int64_t minimumRepresentableExponent = [] {
        const double minimumPositive =
            std::numeric_limits<double>::denorm_min() > 0.0
                ? std::numeric_limits<double>::denorm_min()
                : std::numeric_limits<double>::min();
        return static_cast<std::int64_t>(
            std::ceil(std::log10(minimumPositive)));
    }();
    static const std::int64_t maximumRepresentableExponent =
        static_cast<std::int64_t>(
            std::floor(std::log10(std::numeric_limits<double>::max())));
    if (exponent < minimumRepresentableExponent)
    {
        return 0.0;
    }
    if (exponent > maximumRepresentableExponent)
    {
        return std::numeric_limits<double>::infinity();
    }
    if (exponent >= minimum && exponent <= maximum)
    {
        return values[static_cast<std::size_t>(exponent - minimum)];
    }
    return std::pow(10.0, static_cast<double>(exponent));
}
}

template <typename Function>
RootResult complexSecantImpl(std::complex<double> initialGuess,
                             double tolerance,
                             int maxIterations,
                             const RootConvergenceSpec &spec,
                             const Function &function)
{
    if (!(tolerance > 0.0) || !std::isfinite(tolerance) ||
        !(spec.geometricTolerance > 0.0) ||
        !std::isfinite(spec.geometricTolerance) ||
        !std::isfinite(spec.maximumLog10BackwardError) ||
        !(spec.denominatorRelativeTolerance >= 0.0) ||
        !std::isfinite(spec.denominatorRelativeTolerance))
    {
        return makeFailure(
            initialGuess, initialGuess, 0, RootFailure::InvalidTolerance,
            std::numeric_limits<double>::infinity(),
            std::numeric_limits<double>::infinity());
    }
    if (maxIterations < 1)
    {
        return makeFailure(
            initialGuess, initialGuess, 0, RootFailure::IterationLimit,
            std::numeric_limits<double>::infinity(),
            std::numeric_limits<double>::infinity());
    }

    std::complex<double> x2 = initialGuess;
    std::complex<double> x1 = x2 + 100.0 * tolerance;
    ScaledComplex f1 = function(x1);
    if (!finite(f1.value))
    {
        return makeFailure(
            initialGuess, x2, 0, RootFailure::NonFiniteValue,
            std::numeric_limits<double>::infinity(),
            std::numeric_limits<double>::infinity());
    }

    double lastAbsoluteCorrection =
        std::numeric_limits<double>::infinity();
    ScaledComplex lastSlope;
    bool hasLastSlope = false;
    for (int iteration = 1; iteration <= maxIterations; ++iteration)
    {
        const std::complex<double> x0 = x1;
        const ScaledComplex f0 = f1;
        x1 = x2;
        f1 = function(x1);
        if (!finite(f1.value))
        {
            return makeFailure(
                initialGuess, x1, iteration, RootFailure::NonFiniteValue,
                std::numeric_limits<double>::infinity(),
                lastAbsoluteCorrection);
        }

        const double currentLog10Residual = residualLog10(f1);
        const double priorCorrection = std::abs(x1 - x0);

        // Scale both function values to the larger decimal exponent. The two
        // power10 calls therefore never overflow.
        const int referencePower = std::max(f0.power10, f1.power10);
        const std::int64_t f0PowerOffset =
            static_cast<std::int64_t>(f0.power10) -
            static_cast<std::int64_t>(referencePower);
        const std::int64_t f1PowerOffset =
            static_cast<std::int64_t>(f1.power10) -
            static_cast<std::int64_t>(referencePower);
        const std::complex<double> scaledF0 =
            f0.value * power10(f0PowerOffset);
        const std::complex<double> scaledF1 =
            f1.value * power10(f1PowerOffset);
        const std::complex<double> denominator = scaledF1 - scaledF0;
        const std::complex<double> delta = x1 - x0;
        const std::complex<double> numerator = scaledF1 * delta;
        if (!finite(denominator))
        {
            return makeFailure(
                initialGuess, x1, iteration, RootFailure::NonFiniteValue,
                currentLog10Residual, priorCorrection);
        }
        const bool numeratorOverflow =
            finite(scaledF1) && finite(delta) && !finite(numerator);
        if (!finite(numerator) && !numeratorOverflow)
        {
            return makeFailure(
                initialGuess, x1, iteration, RootFailure::NonFiniteValue,
                currentLog10Residual, priorCorrection);
        }

        const double denominatorScale =
            std::abs(scaledF1) + std::abs(scaledF0);
        const bool exactlyDegenerate =
            denominatorScale == 0.0 || std::abs(denominator) == 0.0;
        const bool relativelyNearDegenerate =
            !exactlyDegenerate &&
            std::abs(denominator) <=
                spec.denominatorRelativeTolerance * denominatorScale;
        if (relativelyNearDegenerate)
        {
            if (!makeScaledSlope(
                    denominator, referencePower, delta, lastSlope))
            {
                return makeFailure(
                    initialGuess, x1, iteration,
                    RootFailure::NonFiniteValue,
                    currentLog10Residual, priorCorrection);
            }
            hasLastSlope = true;
        }
        std::complex<double> shift;
        if (exactlyDegenerate || relativelyNearDegenerate)
        {
            const double currentLog10BackwardError =
                backwardErrorLog10(f1, lastSlope, hasLastSlope);
            const double currentAbsoluteCorrection = std::max(
                priorCorrection,
                magnitudeFromLog10(currentLog10BackwardError));
            const bool geometrySatisfied =
                priorCorrection < spec.geometricTolerance;
            if (geometrySatisfied &&
                currentLog10BackwardError <=
                    spec.maximumLog10BackwardError)
            {
                return makeSuccess(
                    initialGuess, x1, iteration,
                    currentLog10Residual, currentAbsoluteCorrection);
            }
            if (!hasLastSlope && !exactZero(f1))
            {
                return makeFailure(
                    initialGuess, x1, iteration,
                    RootFailure::DegenerateSecant,
                    currentLog10Residual, currentAbsoluteCorrection);
            }
            shift = {0.1 * tolerance, 0.0};
        }
        else
        {
            if (!makeScaledSlope(
                    denominator, referencePower, delta, lastSlope))
            {
                return makeFailure(
                    initialGuess, x1, iteration,
                    RootFailure::NonFiniteValue,
                    currentLog10Residual, priorCorrection);
            }
            hasLastSlope = true;

            // KrakenC bounds any step which is at least as large as x1.
            shift = numeratorOverflow ||
                            std::abs(numerator) >=
                            std::abs(denominator * x1)
                        ? std::complex<double>{0.1 * tolerance, 0.0}
                        : numerator / denominator;
        }

        const std::complex<double> next = x1 - shift;
        if (!finite(next))
        {
            return makeFailure(
                initialGuess, x1, iteration, RootFailure::NonFiniteValue,
                currentLog10Residual, std::abs(shift));
        }

        const double correction = std::abs(next - x1);
        const double threePointCorrection =
            correction + std::abs(next - x0);
        lastAbsoluteCorrection = threePointCorrection;
        const ScaledComplex nextValue = function(next);
        if (!finite(nextValue.value))
        {
            return makeFailure(
                initialGuess, next, iteration, RootFailure::NonFiniteValue,
                std::numeric_limits<double>::infinity(),
                lastAbsoluteCorrection);
        }
        const double nextLog10Residual = residualLog10(nextValue);
        const double nextLog10BackwardError =
            backwardErrorLog10(nextValue, lastSlope, hasLastSlope);
        lastAbsoluteCorrection = std::max(
            threePointCorrection,
            magnitudeFromLog10(nextLog10BackwardError));
        const bool geometrySatisfied =
            threePointCorrection < spec.geometricTolerance;
        if (geometrySatisfied &&
            nextLog10BackwardError <=
                spec.maximumLog10BackwardError)
        {
            return makeSuccess(
                initialGuess, next, iteration,
                nextLog10Residual, lastAbsoluteCorrection);
        }
        x2 = next;
    }

    const ScaledComplex finalValue = function(x2);
    if (finite(finalValue.value))
    {
        lastAbsoluteCorrection = std::max(
            lastAbsoluteCorrection,
            magnitudeFromLog10(
                backwardErrorLog10(finalValue, lastSlope, hasLastSlope)));
    }
    return makeFailure(
        initialGuess, x2, maxIterations, RootFailure::IterationLimit,
        finite(finalValue.value)
            ? residualLog10(finalValue)
            : std::numeric_limits<double>::infinity(),
        lastAbsoluteCorrection);
}

RootResult complexSecant(
    std::complex<double> initialGuess,
    double tolerance,
    int maxIterations,
    const ScaledFunction &function,
    const RootConvergenceSpec &spec)
{
    return complexSecantImpl(
        initialGuess, tolerance, maxIterations, spec, function);
}

RootResult complexSecant(std::complex<double> initialGuess,
                         double tolerance,
                         int maxIterations,
                         const ScaledFunction &function)
{
    return complexSecant(
        initialGuess, tolerance, maxIterations, function,
        defaultConvergenceSpec(tolerance));
}

RootResult complexSecantDispersion(
    std::complex<double> initialGuess,
    double tolerance,
    int maxIterations,
    const ComplexDispersion &dispersion,
    const std::vector<std::complex<double>> &acceptedRoots,
    const RootConvergenceSpec &spec)
{
    const auto function = [&](std::complex<double> value) {
        return dispersion.evaluate(value, acceptedRoots);
    };
    return complexSecantImpl(
        initialGuess, tolerance, maxIterations, spec, function);
}

RootResult complexSecantDispersion(
    std::complex<double> initialGuess,
    double tolerance,
    int maxIterations,
    const ComplexDispersion &dispersion,
    const std::vector<std::complex<double>> &acceptedRoots)
{
    return complexSecantDispersion(
        initialGuess, tolerance, maxIterations,
        dispersion, acceptedRoots,
        defaultConvergenceSpec(tolerance));
}
}
