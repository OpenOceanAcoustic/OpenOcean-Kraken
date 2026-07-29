#include "algorithm/SspInterpolation.h"

#include "algorithm/ComplexNumerics.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>

namespace OpenOceanKrakenc
{
namespace
{
std::string layerError(std::size_t mediumNumber,
                       const std::string &reason)
{
    return "medium " + std::to_string(mediumNumber) + ": " + reason;
}

std::string nodeError(std::size_t mediumNumber,
                      std::size_t sampleNumber,
                      double depth,
                      const char *component,
                      const std::string &reason)
{
    std::ostringstream stream;
    stream << "medium " << mediumNumber
           << " " << component
           << " sample " << sampleNumber
           << " depth " << std::setprecision(17) << depth
           << ": " << reason;
    return stream.str();
}

std::string intervalError(std::size_t mediumNumber,
                          std::size_t intervalNumber,
                          const char *component,
                          const std::string &reason)
{
    return "medium " + std::to_string(mediumNumber) +
           " " + component +
           " interval " + std::to_string(intervalNumber) +
           ": " + reason;
}

std::string nodeRangeError(std::size_t mediumNumber,
                           std::size_t nodeCount,
                           const char *component,
                           const std::string &reason)
{
    return "medium " + std::to_string(mediumNumber) +
           " " + component +
           " nodes 1-" + std::to_string(nodeCount) +
           ": " + reason;
}

bool finiteSample(const AcousticSample &sample)
{
    return std::isfinite(sample.depth) &&
           std::isfinite(sample.cp) &&
           std::isfinite(sample.cs) &&
           std::isfinite(sample.rho) &&
           std::isfinite(sample.alphaP) &&
           std::isfinite(sample.alphaS);
}

bool finiteComplex(std::complex<double> value)
{
    return std::isfinite(value.real()) &&
           std::isfinite(value.imag());
}

const char *unsupportedInterpolationMessage(
    AcousticInterpolation interpolation)
{
    switch (interpolation)
    {
    case AcousticInterpolation::Pchip:
        return "Pchip SSP interpolation is not implemented";
    case AcousticInterpolation::CubicSpline:
        return "cubic-spline SSP interpolation is not implemented";
    case AcousticInterpolation::N2Linear:
    case AcousticInterpolation::CLinear:
        break;
    }
    return "unsupported SSP interpolation";
}

bool supportedInterpolation(AcousticInterpolation interpolation)
{
    return interpolation == AcousticInterpolation::N2Linear ||
           interpolation == AcousticInterpolation::CLinear ||
           interpolation == AcousticInterpolation::CubicSpline ||
           interpolation == AcousticInterpolation::Pchip;
}

std::array<std::complex<double>, 4> linearCoefficients(
    std::complex<double> left,
    std::complex<double> right)
{
    return {left, right - left, std::complex<double>{},
            std::complex<double>{}};
}

std::complex<double> evaluatePolynomial(
    const std::array<std::complex<double>, 4> &coefficients,
    double coordinate)
{
    return coefficients[0] +
           coordinate *
               (coefficients[1] +
                coordinate *
                    (coefficients[2] +
                     coordinate * coefficients[3]));
}

struct SplineRecord
{
    std::complex<double> value{};
    std::complex<double> slope{};
    std::complex<double> curvature{};
    std::complex<double> curvatureRate{};
};

enum class SplineBoundary
{
    NotAKnot,
    FirstDerivative
};

std::vector<std::array<std::complex<double>, 4>> buildSpline(
    const std::vector<double> &coordinates,
    const std::vector<std::complex<double>> &values,
    SplineBoundary leftBoundary,
    std::complex<double> leftSlope,
    SplineBoundary rightBoundary,
    std::complex<double> rightSlope,
    std::size_t mediumNumber,
    const char *component)
{
    const std::size_t count = values.size();
    const auto requireSystemValue =
        [&](std::complex<double> value)
    {
        if (!finiteComplex(value) ||
            value == std::complex<double>{})
        {
            throw std::invalid_argument(nodeRangeError(
                mediumNumber, count, component,
                "singular spline system"));
        }
    };
    std::vector<SplineRecord> work(count);
    for (std::size_t index = 0; index < count; ++index)
    {
        work[index].value = values[index];
    }
    work.front().slope = leftSlope;
    work.back().slope = rightSlope;

    const std::size_t lastInterval = count - 1;
    for (std::size_t index = 1; index < count; ++index)
    {
        work[index].curvature =
            coordinates[index] - coordinates[index - 1];
        requireSystemValue(work[index].curvature);
        work[index].curvatureRate =
            (work[index].value - work[index - 1].value) /
            work[index].curvature;
    }

    if (leftBoundary == SplineBoundary::NotAKnot)
    {
        if (count > 2)
        {
            work[0].curvatureRate = work[2].curvature;
            work[0].curvature =
                work[1].curvature + work[2].curvature;
            requireSystemValue(work[0].curvature);
            work[0].slope =
                ((work[1].curvature +
                  2.0 * work[0].curvature) *
                     work[1].curvatureRate *
                     work[2].curvature +
                 work[1].curvature * work[1].curvature *
                     work[2].curvatureRate) /
                work[0].curvature;
        }
        else
        {
            work[0].curvatureRate = {1.0, 0.0};
            work[0].curvature = {1.0, 0.0};
            work[0].slope =
                2.0 * work[1].curvatureRate;
        }
    }
    else
    {
        work[0].curvatureRate = {1.0, 0.0};
        work[0].curvature = {0.0, 0.0};
    }

    for (std::size_t index = 1;
         index < lastInterval; ++index)
    {
        requireSystemValue(work[index - 1].curvatureRate);
        const std::complex<double> elimination =
            -work[index + 1].curvature /
            work[index - 1].curvatureRate;
        requireSystemValue(elimination);
        work[index].slope =
            elimination * work[index - 1].slope +
            3.0 *
                (work[index].curvature *
                     work[index + 1].curvatureRate +
                 work[index + 1].curvature *
                     work[index].curvatureRate);
        work[index].curvatureRate =
            elimination * work[index - 1].curvature +
            2.0 *
                (work[index].curvature +
                 work[index + 1].curvature);
        requireSystemValue(work[index].curvatureRate);
    }

    if (rightBoundary != SplineBoundary::FirstDerivative)
    {
        std::complex<double> elimination{};
        if (count == 2 &&
            leftBoundary == SplineBoundary::NotAKnot)
        {
            work[count - 1].slope =
                work[count - 1].curvatureRate;
        }
        else if ((count == 3 &&
                  leftBoundary == SplineBoundary::NotAKnot) ||
                 count == 2)
        {
            work[count - 1].slope =
                2.0 * work[count - 1].curvatureRate;
            work[count - 1].curvatureRate =
                {1.0, 0.0};
            requireSystemValue(work[count - 2].curvatureRate);
            elimination =
                -1.0 / work[count - 2].curvatureRate;
            requireSystemValue(elimination);
        }
        else
        {
            elimination =
                work[count - 2].curvature +
                work[count - 1].curvature;
            requireSystemValue(elimination);
            work[count - 1].slope =
                ((work[count - 1].curvature +
                  2.0 * elimination) *
                     work[count - 1].curvatureRate *
                     work[count - 2].curvature +
                 work[count - 1].curvature *
                     work[count - 1].curvature *
                     (work[count - 2].value -
                      work[count - 3].value) /
                     work[count - 2].curvature) /
                elimination;
            requireSystemValue(work[count - 2].curvatureRate);
            elimination =
                -elimination /
                work[count - 2].curvatureRate;
            requireSystemValue(elimination);
            work[count - 1].curvatureRate =
                work[count - 2].curvature;
        }

        if (leftBoundary ==
                SplineBoundary::FirstDerivative ||
            count > 2)
        {
            work[count - 1].curvatureRate =
                elimination * work[count - 2].curvature +
                work[count - 1].curvatureRate;
            requireSystemValue(work[count - 1].curvatureRate);
            work[count - 1].slope =
                (elimination * work[count - 2].slope +
                 work[count - 1].slope) /
                work[count - 1].curvatureRate;
        }
    }

    for (std::size_t index = lastInterval;
         index-- > 0;)
    {
        requireSystemValue(work[index].curvatureRate);
        work[index].slope =
            (work[index].slope -
             work[index].curvature *
                 work[index + 1].slope) /
            work[index].curvatureRate;
    }

    for (std::size_t index = 1; index < count; ++index)
    {
        const std::complex<double> delta =
            work[index].curvature;
        const std::complex<double> firstDifference =
            (work[index].value -
             work[index - 1].value) /
            delta;
        const std::complex<double> thirdDifference =
            work[index - 1].slope +
            work[index].slope -
            2.0 * firstDifference;
        work[index - 1].curvature =
            2.0 *
            (firstDifference -
             work[index - 1].slope -
             thirdDifference) /
            delta;
        work[index - 1].curvatureRate =
            (thirdDifference / delta) *
            (6.0 / delta);
    }

    work[count - 1].curvature =
        work[lastInterval - 1].curvature +
        (coordinates[count - 1] -
         coordinates[lastInterval - 1]) *
            work[lastInterval - 1].curvatureRate;

    work[count - 1].curvatureRate = {};
    for (std::size_t index = 0;
         index < lastInterval; ++index)
    {
        const std::complex<double> delta =
            coordinates[index + 1] - coordinates[index];
        work[count - 1].curvatureRate +=
            delta *
            (work[index].value +
             delta *
                 (work[index].slope / 2.0 +
                  delta *
                      (work[index].curvature / 6.0 +
                       delta *
                           work[index].curvatureRate /
                           24.0)));
    }
    work[count - 1].curvatureRate /=
        coordinates[count - 1] - coordinates[0];

    std::vector<std::array<std::complex<double>, 4>>
        segments;
    segments.reserve(lastInterval);
    for (std::size_t index = 0;
         index < lastInterval; ++index)
    {
        segments.push_back(
            {work[index].value,
             work[index].slope,
             work[index].curvature / 2.0,
             work[index].curvatureRate / 6.0});
    }
    return segments;
}

double pchipLeftEndpointDerivative(double delta1,
                                   double delta2,
                                   double derivative)
{
    if (delta1 * derivative <= 0.0)
    {
        return 0.0;
    }
    if (delta1 * delta2 <= 0.0 &&
        std::abs(derivative) > std::abs(3.0 * delta1))
    {
        return 3.0 * delta1;
    }
    return derivative;
}

double pchipRightEndpointDerivative(double delta1,
                                    double delta2,
                                    double derivative)
{
    if (delta2 * derivative <= 0.0)
    {
        return 0.0;
    }
    if (delta1 * delta2 <= 0.0 &&
        std::abs(derivative) > std::abs(3.0 * delta2))
    {
        return 3.0 * delta2;
    }
    return derivative;
}

double pchipInteriorDerivative(double delta1,
                               double delta2,
                               double derivative)
{
    if (delta1 * delta2 <= 0.0)
    {
        return 0.0;
    }
    if (delta1 > 0.0)
    {
        return std::min(
            std::max(derivative, 0.0),
            3.0 * std::min(delta1, delta2));
    }
    return std::max(
        std::min(derivative, 0.0),
        3.0 * std::max(delta1, delta2));
}

std::complex<double> pchipLeftEndpointDerivative(
    std::complex<double> delta1,
    std::complex<double> delta2,
    std::complex<double> derivative)
{
    return {
        pchipLeftEndpointDerivative(
            delta1.real(), delta2.real(), derivative.real()),
        pchipLeftEndpointDerivative(
            delta1.imag(), delta2.imag(), derivative.imag())};
}

std::complex<double> pchipRightEndpointDerivative(
    std::complex<double> delta1,
    std::complex<double> delta2,
    std::complex<double> derivative)
{
    return {
        pchipRightEndpointDerivative(
            delta1.real(), delta2.real(), derivative.real()),
        pchipRightEndpointDerivative(
            delta1.imag(), delta2.imag(), derivative.imag())};
}

std::complex<double> pchipInteriorDerivative(
    std::complex<double> delta1,
    std::complex<double> delta2,
    std::complex<double> derivative)
{
    return {
        pchipInteriorDerivative(
            delta1.real(), delta2.real(), derivative.real()),
        pchipInteriorDerivative(
            delta1.imag(), delta2.imag(), derivative.imag())};
}

std::vector<std::array<std::complex<double>, 4>> buildPchip(
    const std::vector<double> &coordinates,
    const std::vector<std::complex<double>> &values,
    std::size_t mediumNumber,
    const char *component)
{
    const std::size_t count = values.size();
    const std::size_t segmentCount = count - 1;
    if (count == 2)
    {
        const double interval =
            coordinates[1] - coordinates[0];
        return {{
            values[0],
            (values[1] - values[0]) / interval,
            std::complex<double>{},
            std::complex<double>{}}};
    }

    std::vector<std::complex<double>> derivatives(count);
    const double leftH1 = coordinates[1] - coordinates[0];
    const double leftH2 = coordinates[2] - coordinates[1];
    const std::complex<double> leftDelta1 =
        (values[1] - values[0]) / leftH1;
    const std::complex<double> leftDelta2 =
        (values[2] - values[1]) / leftH2;
    const std::complex<double> leftRaw =
        ((2.0 * leftH1 + leftH2) * leftDelta1 -
         leftH1 * leftDelta2) /
        (leftH1 + leftH2);
    derivatives.front() = pchipLeftEndpointDerivative(
        leftDelta1, leftDelta2, leftRaw);

    const double rightH1 =
        coordinates[count - 2] - coordinates[count - 3];
    const double rightH2 =
        coordinates[count - 1] - coordinates[count - 2];
    const std::complex<double> rightDelta1 =
        (values[count - 2] - values[count - 3]) / rightH1;
    const std::complex<double> rightDelta2 =
        (values[count - 1] - values[count - 2]) / rightH2;
    const std::complex<double> rightRaw =
        (-rightH2 * rightDelta1 +
         (rightH1 + 2.0 * rightH2) * rightDelta2) /
        (rightH1 + rightH2);
    derivatives.back() = pchipRightEndpointDerivative(
        rightDelta1, rightDelta2, rightRaw);

    const auto spline = buildSpline(
        coordinates, values,
        SplineBoundary::FirstDerivative, derivatives.front(),
        SplineBoundary::FirstDerivative, derivatives.back(),
        mediumNumber, component);
    for (std::size_t index = 1;
         index + 1 < count; ++index)
    {
        const double leftInterval =
            coordinates[index] - coordinates[index - 1];
        const double rightInterval =
            coordinates[index + 1] - coordinates[index];
        const std::complex<double> delta1 =
            (values[index] - values[index - 1]) /
            leftInterval;
        const std::complex<double> delta2 =
            (values[index + 1] - values[index]) /
            rightInterval;
        derivatives[index] = pchipInteriorDerivative(
            delta1, delta2, spline[index][1]);
    }

    std::vector<std::array<std::complex<double>, 4>>
        segments;
    segments.reserve(segmentCount);
    for (std::size_t index = 0;
         index < segmentCount; ++index)
    {
        const double interval =
            coordinates[index + 1] - coordinates[index];
        const std::complex<double> difference =
            values[index + 1] - values[index];
        segments.push_back({
            values[index],
            derivatives[index],
            (3.0 * difference -
             interval *
                 (2.0 * derivatives[index] +
                  derivatives[index + 1])) /
                (interval * interval),
            (interval *
                 (derivatives[index] +
                  derivatives[index + 1]) -
             2.0 * difference) /
                (interval * interval * interval)});
    }
    return segments;
}
}

void validateSspLayer(const AcousticLayer &layer,
                      std::size_t mediumNumber)
{
    if (layer.baseMesh < 1)
    {
        throw std::invalid_argument(
            layerError(mediumNumber, "base mesh must be positive"));
    }
    if (!std::isfinite(layer.topDepth) ||
        !std::isfinite(layer.bottomDepth) ||
        !std::isfinite(layer.roughnessRms) ||
        !std::isfinite(layer.attenuationPower) ||
        !std::isfinite(layer.transitionFrequency))
    {
        throw std::invalid_argument(layerError(
            mediumNumber,
            "layer bounds, roughness, attenuation power, and transition frequency must be finite"));
    }
    if (!(layer.bottomDepth > layer.topDepth))
    {
        throw std::invalid_argument(layerError(
            mediumNumber,
            "layer bottom depth must be greater than top depth"));
    }
    if (layer.roughnessRms < 0.0)
    {
        throw std::invalid_argument(
            layerError(mediumNumber,
                       "layer roughness must be non-negative"));
    }

    if (layer.samples.size() < 2)
    {
        throw std::invalid_argument(layerError(
            mediumNumber,
            "acoustic layer requires at least two SSP samples"));
    }

    constexpr double endpointTolerance =
        100.0 * std::numeric_limits<float>::epsilon();
    if (!(std::abs(layer.samples.front().depth -
                   layer.topDepth) < endpointTolerance) ||
        !(std::abs(layer.samples.back().depth -
                   layer.bottomDepth) < endpointTolerance))
    {
        throw std::invalid_argument(layerError(
            mediumNumber,
            "SSP samples must span the layer bounds"));
    }

    for (std::size_t index = 0;
         index < layer.samples.size(); ++index)
    {
        const AcousticSample &sample = layer.samples[index];
        if (!finiteSample(sample))
        {
            throw std::invalid_argument(nodeError(
                mediumNumber, index + 1, sample.depth,
                "sample", "all fields must be finite"));
        }
        if (index > 0 &&
            !(sample.depth >
              layer.samples[index - 1].depth))
        {
            throw std::invalid_argument(nodeError(
                mediumNumber, index + 1, sample.depth,
                "depth", "SSP depths must be strictly increasing"));
        }
        if (!(sample.cp > 0.0))
        {
            throw std::invalid_argument(nodeError(
                mediumNumber, index + 1, sample.depth,
                "cp", "sound speed must be positive"));
        }
        if (!(sample.rho > 0.0))
        {
            throw std::invalid_argument(nodeError(
                mediumNumber, index + 1, sample.depth,
                "rho", "density must be positive"));
        }
        if (sample.alphaP < 0.0)
        {
            throw std::invalid_argument(nodeError(
                mediumNumber, index + 1, sample.depth,
                "cp", "attenuation must be non-negative"));
        }
        if (sample.alphaS < 0.0)
        {
            throw std::invalid_argument(nodeError(
                mediumNumber, index + 1, sample.depth,
                "cs", "attenuation must be non-negative"));
        }
    }

    const bool elastic = layer.samples.front().cs > 0.0;
    for (std::size_t index = 0;
         index < layer.samples.size(); ++index)
    {
        const AcousticSample &sample = layer.samples[index];
        if (!(elastic ? sample.cs > 0.0
                      : sample.cs == 0.0))
        {
            throw std::invalid_argument(nodeError(
                mediumNumber, index + 1, sample.depth, "cs",
                "sound speeds must be either all zero or all positive"));
        }
    }
    if (elastic && layer.roughnessRms != 0.0)
    {
        throw std::invalid_argument(
            "Rough elastic interfaces are not allowed");
    }
}

PreparedComplexSspLayer::PreparedComplexSspLayer(
    const AcousticLayer &layer,
    AcousticInterpolation interpolation,
    const AttenuationContext &context,
    std::size_t mediumNumber)
    : interpolation_(interpolation),
      mediumNumber_(mediumNumber)
{
    validateSspLayer(layer, mediumNumber_);
    if (!supportedInterpolation(interpolation_))
    {
        throw std::invalid_argument(
            unsupportedInterpolationMessage(interpolation_));
    }

    elastic_ = layer.samples.front().cs > 0.0;
    depths_.reserve(layer.samples.size());
    nodes_.reserve(layer.samples.size());
    for (std::size_t index = 0;
         index < layer.samples.size(); ++index)
    {
        const AcousticSample &sample = layer.samples[index];
        const auto complexify =
            [&](const char *component,
                double speed,
                double attenuation)
        {
            try
            {
                return complexSoundSpeed(
                    sample.depth, speed, attenuation, context);
            }
            catch (const std::exception &error)
            {
                throw std::invalid_argument(nodeError(
                    mediumNumber_, index + 1, sample.depth,
                    component, error.what()));
            }
        };

        InterpolatedMaterial node;
        node.cp = complexify("cp", sample.cp, sample.alphaP);
        node.cs = elastic_
                      ? complexify("cs", sample.cs, sample.alphaS)
                      : std::complex<double>{};
        node.rho = sample.rho;
        depths_.push_back(sample.depth);
        nodes_.push_back(node);
    }

    const std::size_t segmentCount = nodes_.size() - 1;
    if (interpolation_ ==
            AcousticInterpolation::CubicSpline ||
        interpolation_ == AcousticInterpolation::Pchip)
    {
        std::vector<std::complex<double>> cpValues;
        std::vector<std::complex<double>> csValues;
        std::vector<std::complex<double>> rhoValues;
        cpValues.reserve(nodes_.size());
        csValues.reserve(nodes_.size());
        rhoValues.reserve(nodes_.size());
        for (const InterpolatedMaterial &node : nodes_)
        {
            cpValues.push_back(node.cp);
            if (elastic_)
            {
                csValues.push_back(node.cs);
            }
            rhoValues.emplace_back(node.rho, 0.0);
        }

        cpCoefficients_ =
            interpolation_ == AcousticInterpolation::Pchip
                ? buildPchip(
                      depths_, cpValues, mediumNumber_, "cp")
                : buildSpline(
                      depths_, cpValues,
                      SplineBoundary::NotAKnot, {},
                      SplineBoundary::NotAKnot, {},
                      mediumNumber_, "cp");
        csCoefficients_ =
            elastic_
                ? (interpolation_ == AcousticInterpolation::Pchip
                       ? buildPchip(
                             depths_, csValues,
                             mediumNumber_, "cs")
                       : buildSpline(
                             depths_, csValues,
                             SplineBoundary::NotAKnot, {},
                             SplineBoundary::NotAKnot, {},
                             mediumNumber_, "cs"))
                : std::vector<SegmentCoefficients>(
                      segmentCount);
        rhoCoefficients_ =
            interpolation_ == AcousticInterpolation::Pchip
                ? buildPchip(
                      depths_, rhoValues, mediumNumber_, "rho")
                : buildSpline(
                      depths_, rhoValues,
                      SplineBoundary::NotAKnot, {},
                      SplineBoundary::NotAKnot, {},
                      mediumNumber_, "rho");
    }
    else
    {
        cpCoefficients_.reserve(segmentCount);
        csCoefficients_.reserve(segmentCount);
        rhoCoefficients_.reserve(segmentCount);
        for (std::size_t index = 0;
             index < segmentCount; ++index)
        {
            if (interpolation_ ==
                AcousticInterpolation::CLinear)
            {
                cpCoefficients_.push_back(
                    linearCoefficients(
                        nodes_[index].cp,
                        nodes_[index + 1].cp));
                csCoefficients_.push_back(
                    elastic_
                        ? linearCoefficients(
                              nodes_[index].cs,
                              nodes_[index + 1].cs)
                        : SegmentCoefficients{});
            }
            else
            {
                const std::complex<double> leftCp =
                    1.0 / (nodes_[index].cp *
                           nodes_[index].cp);
                const std::complex<double> rightCp =
                    1.0 / (nodes_[index + 1].cp *
                           nodes_[index + 1].cp);
                cpCoefficients_.push_back(
                    linearCoefficients(leftCp, rightCp));
                if (elastic_)
                {
                    const std::complex<double> leftCs =
                        1.0 / (nodes_[index].cs *
                               nodes_[index].cs);
                    const std::complex<double> rightCs =
                        1.0 /
                        (nodes_[index + 1].cs *
                         nodes_[index + 1].cs);
                    csCoefficients_.push_back(
                        linearCoefficients(
                            leftCs, rightCs));
                }
                else
                {
                    csCoefficients_.push_back(
                        SegmentCoefficients{});
                }
            }
            rhoCoefficients_.push_back(
                linearCoefficients(
                    nodes_[index].rho,
                    nodes_[index + 1].rho));
        }
    }

    const auto validateCoefficients =
        [&](const std::vector<SegmentCoefficients> &segments,
            const char *component)
    {
        for (std::size_t index = 0;
             index < segments.size(); ++index)
        {
            for (const std::complex<double> coefficient :
                 segments[index])
            {
                if (!finiteComplex(coefficient))
                {
                    throw std::invalid_argument(intervalError(
                        mediumNumber_, index + 1, component,
                        "non-finite spline coefficient"));
                }
            }
        }
    };
    validateCoefficients(cpCoefficients_, "cp");
    validateCoefficients(csCoefficients_, "cs");
    validateCoefficients(rhoCoefficients_, "rho");
}

InterpolatedMaterial PreparedComplexSspLayer::evaluate(
    double depth) const
{
    if (!std::isfinite(depth) ||
        depth < depths_.front() ||
        depth > depths_.back())
    {
        throw std::invalid_argument(layerError(
            mediumNumber_,
            "SSP query depth is outside the layer"));
    }

    const auto upper = std::lower_bound(
        depths_.begin(), depths_.end(), depth);
    const std::size_t upperIndex =
        static_cast<std::size_t>(upper - depths_.begin());
    if (upper != depths_.end() && *upper == depth)
    {
        return nodes_[upperIndex];
    }

    const std::size_t leftIndex = upperIndex - 1;
    const double fraction =
        (depth - depths_[leftIndex]) /
        (depths_[leftIndex + 1] - depths_[leftIndex]);
    const double coordinate =
        (interpolation_ == AcousticInterpolation::CubicSpline ||
         interpolation_ == AcousticInterpolation::Pchip)
            ? depth - depths_[leftIndex]
            : fraction;
    const std::complex<double> cpInterpolant =
        evaluatePolynomial(
            cpCoefficients_[leftIndex], coordinate);
    const std::complex<double> csInterpolant =
        evaluatePolynomial(
            csCoefficients_[leftIndex], coordinate);
    const std::complex<double> rhoInterpolant =
        evaluatePolynomial(
            rhoCoefficients_[leftIndex], coordinate);

    InterpolatedMaterial result;
    if (interpolation_ == AcousticInterpolation::CLinear)
    {
        result.cp =
            (1.0 - fraction) * nodes_[leftIndex].cp +
            fraction * nodes_[leftIndex + 1].cp;
        result.cs =
            elastic_
                ? (1.0 - fraction) * nodes_[leftIndex].cs +
                      fraction * nodes_[leftIndex + 1].cs
                : std::complex<double>{};
    }
    else if (interpolation_ ==
                 AcousticInterpolation::CubicSpline ||
        interpolation_ == AcousticInterpolation::Pchip)
    {
        result.cp = cpInterpolant;
        result.cs =
            elastic_ ? csInterpolant
                     : std::complex<double>{};
    }
    else if (interpolation_ ==
             AcousticInterpolation::N2Linear)
    {
        result.cp =
            1.0 / std::sqrt(cpInterpolant);
        result.cs =
            elastic_
                ? 1.0 / std::sqrt(csInterpolant)
                : std::complex<double>{};
    }
    else
    {
        throw std::invalid_argument(
            unsupportedInterpolationMessage(interpolation_));
    }
    if ((interpolation_ ==
             AcousticInterpolation::CubicSpline ||
         interpolation_ == AcousticInterpolation::Pchip) &&
        std::abs(rhoInterpolant.imag()) > 1.0e-13)
    {
        throw std::invalid_argument(layerError(
            mediumNumber_,
            "density spline produced a non-real value"));
    }
    result.rho =
        interpolation_ == AcousticInterpolation::CLinear
            ? (1.0 - fraction) * nodes_[leftIndex].rho +
                  fraction * nodes_[leftIndex + 1].rho
            : rhoInterpolant.real();
    return result;
}

const InterpolatedMaterial &
PreparedComplexSspLayer::top() const noexcept
{
    return nodes_.front();
}

const InterpolatedMaterial &
PreparedComplexSspLayer::bottom() const noexcept
{
    return nodes_.back();
}

bool PreparedComplexSspLayer::elastic() const noexcept
{
    return elastic_;
}

double interpolateRealCp(const AcousticLayer &layer,
                         AcousticInterpolation interpolation,
                         double depth,
                         std::size_t mediumNumber)
{
    validateSspLayer(layer, mediumNumber);
    if (!supportedInterpolation(interpolation))
    {
        throw std::invalid_argument(
            unsupportedInterpolationMessage(interpolation));
    }
    if (!std::isfinite(depth) ||
        depth < layer.samples.front().depth ||
        depth > layer.samples.back().depth)
    {
        throw std::invalid_argument(layerError(
            mediumNumber,
            "SSP query depth is outside the layer"));
    }

    const auto upper = std::lower_bound(
        layer.samples.begin(), layer.samples.end(), depth,
        [](const AcousticSample &sample, double query)
        {
            return sample.depth < query;
        });
    if (upper != layer.samples.end() &&
        upper->depth == depth)
    {
        return upper->cp;
    }

    const AcousticSample &right = *upper;
    const AcousticSample &left = *(upper - 1);
    const double fraction =
        (depth - left.depth) /
        (right.depth - left.depth);
    if (interpolation == AcousticInterpolation::CLinear)
    {
        return (1.0 - fraction) * left.cp +
               fraction * right.cp;
    }
    if (interpolation == AcousticInterpolation::N2Linear)
    {
        const double inverseSquared =
            (1.0 - fraction) / (left.cp * left.cp) +
            fraction / (right.cp * right.cp);
        return 1.0 / std::sqrt(inverseSquared);
    }

    std::vector<double> depths;
    std::vector<std::complex<double>> values;
    depths.reserve(layer.samples.size());
    values.reserve(layer.samples.size());
    for (const AcousticSample &sample : layer.samples)
    {
        depths.push_back(sample.depth);
        values.emplace_back(sample.cp, 0.0);
    }
    const auto coefficients =
        interpolation == AcousticInterpolation::Pchip
            ? buildPchip(
                  depths, values, mediumNumber, "cp")
            : buildSpline(
                  depths, values,
                  SplineBoundary::NotAKnot, {},
                  SplineBoundary::NotAKnot, {},
                  mediumNumber, "cp");
    const std::size_t leftIndex =
        static_cast<std::size_t>(
            upper - layer.samples.begin() - 1);
    const std::complex<double> interpolated =
        evaluatePolynomial(
            coefficients[leftIndex],
            depth - left.depth);
    if (std::abs(interpolated.imag()) > 1.0e-13)
    {
        throw std::invalid_argument(layerError(
            mediumNumber,
            "real cp spline produced a non-real value"));
    }
    return interpolated.real();
}
}
