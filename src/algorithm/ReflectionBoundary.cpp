#include "algorithm/ReflectionBoundary.h"

#include "algorithm/ComplexNumerics.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace OpenOceanKrakenc
{
namespace
{
std::complex<double> polynomial(
    std::complex<double> x,
    const std::vector<double> &abscissas,
    const std::vector<std::complex<double>> &values)
{
    std::complex<double> result = 0.0;
    for (std::size_t i = 0; i < values.size(); ++i)
    {
        std::complex<double> term = values[i];
        for (std::size_t j = 0; j < values.size(); ++j)
        {
            if (i != j)
            {
                term *= (x - abscissas[j]) / (abscissas[i] - abscissas[j]);
            }
        }
        result += term;
    }
    return result;
}

Impedance internalImpedance(const AcousticBoundary &boundary,
                            std::complex<double> eigenvalue)
{
    const auto &samples = boundary.internalSamples;
    if (samples.size() < 2)
    {
        throw std::invalid_argument("IRC boundary requires at least two samples");
    }
    if (eigenvalue.real() <= samples.front().eigenvalue)
    {
        return {samples.front().f, samples.front().g, samples.front().power10};
    }
    if (eigenvalue.real() >= samples.back().eigenvalue)
    {
        return {samples.back().f, samples.back().g, samples.back().power10};
    }
    const auto right = std::upper_bound(
        samples.begin(), samples.end(), eigenvalue.real(),
        [](double value, const InternalReflectionSample &sample) {
            return value < sample.eigenvalue;
        });
    std::size_t rightIndex = static_cast<std::size_t>(right - samples.begin());
    std::size_t leftIndex = rightIndex - 1;
    if (leftIndex > 0)
    {
        --leftIndex;
    }
    rightIndex = std::min(leftIndex + 2, samples.size() - 1);
    if (rightIndex - leftIndex + 1 < 3 && leftIndex > 0)
    {
        --leftIndex;
    }

    const int basePower = samples[leftIndex].power10;
    std::vector<double> x;
    std::vector<std::complex<double>> f;
    std::vector<std::complex<double>> g;
    for (std::size_t index = leftIndex; index <= rightIndex; ++index)
    {
        const double scale = std::pow(
            10.0, samples[index].power10 - basePower);
        x.push_back(samples[index].eigenvalue);
        f.push_back(samples[index].f * scale);
        g.push_back(samples[index].g * scale);
    }
    return {polynomial(eigenvalue, x, f), polynomial(eigenvalue, x, g), basePower};
}
}

Impedance reflectionBoundaryImpedance(
    const AcousticBoundary &boundary,
    bool top,
    double omega2,
    std::complex<double> eigenvalue,
    std::complex<double> insideSoundSpeed)
{
    Impedance result;
    if (boundary.type == AcousticBoundaryType::InternalReflection)
    {
        result = internalImpedance(boundary, eigenvalue);
    }
    else if (boundary.type == AcousticBoundaryType::ReflectionCoefficient)
    {
        if (boundary.reflectionSamples.size() < 2)
        {
            throw std::invalid_argument("BRC boundary requires at least two samples");
        }
        const std::complex<double> kx = std::sqrt(eigenvalue);
        const std::complex<double> kz = std::sqrt(
            omega2 / (insideSoundSpeed * insideSoundSpeed) - eigenvalue);
        const double angle = std::atan2(kz.real(), kx.real()) * 180.0 / pi;
        double magnitude = 0.0;
        double phase = 0.0;
        const auto &samples = boundary.reflectionSamples;
        if (angle >= samples.front().angleDegrees &&
            angle <= samples.back().angleDegrees)
        {
            const auto right = std::upper_bound(
                samples.begin(), samples.end(), angle,
                [](double value, const ReflectionSample &sample) {
                    return value < sample.angleDegrees;
                });
            if (right == samples.end())
            {
                magnitude = samples.back().magnitude;
                phase = samples.back().phaseRadians;
            }
            else
            {
                const ReflectionSample &left = *(right - 1);
                const double fraction = (angle - left.angleDegrees) /
                                        (right->angleDegrees - left.angleDegrees);
                magnitude = (1.0 - fraction) * left.magnitude +
                            fraction * right->magnitude;
                phase = (1.0 - fraction) * left.phaseRadians +
                        fraction * right->phaseRadians;
            }
        }
        const std::complex<double> reflection = magnitude *
                                                std::exp(std::complex<double>(0.0, phase));
        if (std::abs(1.0 - reflection) <= 1.0e-13)
        {
            result.f = 0.0;
            result.g = 1.0;
        }
        else
        {
            result.f = 1.0;
            result.g = (1.0 + reflection) /
                       (std::complex<double>(0.0, 1.0) * kz *
                        (1.0 - reflection));
        }
    }
    else
    {
        throw std::invalid_argument("boundary is not a BRC or IRC boundary");
    }
    if (top)
    {
        result.g = -result.g;
    }
    return result;
}
}
