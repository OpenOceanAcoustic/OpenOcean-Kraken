#include "algorithm/ComplexModeNormalization.h"

#include "algorithm/ComplexDispersion.h"
#include "algorithm/ComplexNumerics.h"
#include "algorithm/ElasticCompound.h"
#include "algorithm/ScatteringLoss.h"

#include <cmath>
#include <stdexcept>

namespace OpenOceanKrakenc
{
namespace
{
std::complex<double> boundaryAdmittanceDerivative(
    const AcousticBoundary &boundary,
    bool top,
    const AcousticCase &input,
    const AcousticMatrix &matrix,
    std::complex<double> eigenvalue)
{
    if (boundary.type == AcousticBoundaryType::Vacuum ||
        boundary.type == AcousticBoundaryType::Rigid)
    {
        return 0.0;
    }
    const std::complex<double> x1 = 0.9999999 * eigenvalue;
    const std::complex<double> x2 = 1.0000001 * eigenvalue;
    const Impedance first = caseBoundaryImpedance(input, matrix, top, x1);
    if (first.g == std::complex<double>(0.0, 0.0))
    {
        return 0.0;
    }
    const Impedance second = caseBoundaryImpedance(input, matrix, top, x2);
    return (second.f / second.g - first.f / first.g) / (x2 - x1);
}
}

std::complex<double> interfacialScatterPerturbation(
    const AcousticCase &input,
    const AcousticMatrix &matrix,
    std::complex<double> eigenvalue,
    const Eigen::Ref<const Eigen::VectorXcd> &normalizedMode)
{
    std::size_t firstAcoustic = input.layers.size();
    std::size_t lastAcoustic = input.layers.size();
    for (std::size_t medium = 0; medium < input.layers.size(); ++medium)
    {
        if (!matrix.layerElastic[medium])
        {
            if (firstAcoustic == input.layers.size()) firstAcoustic = medium;
            lastAcoustic = medium;
        }
    }
    if (firstAcoustic == input.layers.size()) return 0.0;

    std::complex<double> perturbation = 0.0;
    int pressureIndex = 0;
    const auto acousticEta = [&](std::size_t medium, int local) {
        const int point = matrix.offsets[medium] + local;
        const double spacing = matrix.spacing[medium];
        return (2.0 + matrix.b1[point]) / (spacing * spacing) - eigenvalue;
    };

    const double topSigma = input.layers[firstAcoustic].roughnessRms;
    if (topSigma != 0.0)
    {
        double rho1 = 0.0;
        std::complex<double> eta1Squared = 0.0;
        std::complex<double> u = 0.0;
        switch (input.top.type)
        {
        case AcousticBoundaryType::HalfSpace:
        {
            const std::complex<double> cp = complexSoundSpeed(
                input.top.depth, input.top.cp, input.top.alphaP,
                attenuationContext(input, input.top.attenuationPower,
                                   input.top.transitionFrequency));
            rho1 = input.top.rho;
            eta1Squared = eigenvalue - matrix.omega2 / (cp * cp);
            u = pekerisRoot(eta1Squared) * normalizedMode[0] / rho1;
            break;
        }
        case AcousticBoundaryType::Vacuum:
            rho1 = 1.0e-9;
            eta1Squared = 1.0;
            u = normalizedMode.size() > 1
                    ? normalizedMode[1] /
                          (matrix.spacing[firstAcoustic] *
                           matrix.rho[matrix.offsets[firstAcoustic]])
                    : std::complex<double>{};
            break;
        case AcousticBoundaryType::Rigid:
            rho1 = 1.0e9;
            eta1Squared = 1.0;
            u = 0.0;
            break;
        default:
            throw std::invalid_argument(
                "rough tabulated top boundaries are not supported by Kuperman-Ingenito");
        }
        const double rho2 = matrix.rho[matrix.offsets[firstAcoustic]];
        perturbation += kupermanIngenito(
            topSigma, eta1Squared, rho1, acousticEta(firstAcoustic, 0), rho2,
            normalizedMode[0], u);
    }

    for (std::size_t medium = firstAcoustic; medium < lastAcoustic; ++medium)
    {
        pressureIndex += matrix.counts[medium];
        const double sigma = input.layers[medium + 1].roughnessRms;
        if (sigma == 0.0) continue;
        if (matrix.layerElastic[medium + 1])
        {
            throw std::invalid_argument("Rough elastic interfaces are not allowed");
        }
        const int upperPoint = matrix.offsets[medium] + matrix.counts[medium];
        const double upperSpacing = matrix.spacing[medium];
        const double upperSpacingSquared = upperSpacing * upperSpacing;
        const double rho1 = matrix.rho[upperPoint];
        const std::complex<double> u =
            (-normalizedMode[pressureIndex - 1] -
             0.5 * (matrix.b1[upperPoint] -
                    upperSpacingSquared * eigenvalue) *
                 normalizedMode[pressureIndex]) /
            (upperSpacing * rho1);
        perturbation += kupermanIngenito(
            sigma, acousticEta(medium, matrix.counts[medium]), rho1,
            acousticEta(medium + 1, 0),
            matrix.rho[matrix.offsets[medium + 1]],
            normalizedMode[pressureIndex], u);
    }

    pressureIndex += matrix.counts[lastAcoustic];
    const double bottomSigma = input.bottom.roughnessRms;
    if (bottomSigma != 0.0)
    {
        const int upperPoint = matrix.offsets[lastAcoustic] +
                               matrix.counts[lastAcoustic];
        const double spacing = matrix.spacing[lastAcoustic];
        const double spacingSquared = spacing * spacing;
        const double rho1 = matrix.rho[upperPoint];
        const std::complex<double> u =
            (-normalizedMode[pressureIndex - 1] -
             0.5 * (matrix.b1[upperPoint] - spacingSquared * eigenvalue) *
                 normalizedMode[pressureIndex]) /
            (spacing * rho1);
        double rho2 = 0.0;
        std::complex<double> eta2Squared = 0.0;
        switch (input.bottom.type)
        {
        case AcousticBoundaryType::HalfSpace:
        {
            const std::complex<double> cp = complexSoundSpeed(
                input.bottom.depth, input.bottom.cp, input.bottom.alphaP,
                attenuationContext(input, input.bottom.attenuationPower,
                                   input.bottom.transitionFrequency));
            rho2 = input.bottom.rho;
            eta2Squared = matrix.omega2 / (cp * cp) - eigenvalue;
            break;
        }
        case AcousticBoundaryType::Vacuum:
            rho2 = 1.0e-9;
            eta2Squared = 1.0;
            break;
        case AcousticBoundaryType::Rigid:
            rho2 = 1.0e9;
            eta2Squared = 1.0;
            break;
        default:
            throw std::invalid_argument(
                "rough tabulated bottom boundaries are not supported by Kuperman-Ingenito");
        }
        perturbation += kupermanIngenito(
            bottomSigma,
            acousticEta(lastAcoustic, matrix.counts[lastAcoustic]), rho1,
            eta2Squared, rho2, normalizedMode[pressureIndex], u);
    }
    return perturbation;
}

namespace
{

std::complex<double> halfSpaceSlowness(
    const AcousticBoundary &boundary,
    const AcousticCase &input,
    const AcousticMatrix &matrix,
    std::complex<double> eigenvalue,
    std::complex<double> pressure)
{
    if (boundary.type != AcousticBoundaryType::HalfSpace)
    {
        return 0.0;
    }
    const std::complex<double> cp = complexSoundSpeed(
        boundary.depth, boundary.cp, boundary.alphaP,
        attenuationContext(input, boundary.attenuationPower,
                           boundary.transitionFrequency));
    const std::complex<double> gamma =
        pekerisRoot(eigenvalue - matrix.omega2 / (cp * cp));
    return pressure * pressure /
           (2.0 * gamma * boundary.rho * cp * cp);
}
}

NormalizedComplexMode normalizeAcousticMode(
    const AcousticCase &input,
    const AcousticMatrix &matrix,
    std::complex<double> eigenvalue,
    int turningPoint,
    const Eigen::Ref<const Eigen::VectorXcd> &rawMode)
{
    std::size_t firstAcoustic = input.layers.size();
    std::size_t lastAcoustic = input.layers.size();
    for (std::size_t medium = 0; medium < input.layers.size(); ++medium)
    {
        if (!matrix.layerElastic[medium])
        {
            if (firstAcoustic == input.layers.size())
            {
                firstAcoustic = medium;
            }
            lastAcoustic = medium;
        }
    }
    int expectedSize = 1;
    for (std::size_t medium = firstAcoustic; medium <= lastAcoustic; ++medium)
    {
        expectedSize += matrix.counts[medium];
    }
    if (rawMode.size() != expectedSize || turningPoint < 0 ||
        turningPoint >= expectedSize || matrix.counts.size() != input.layers.size())
    {
        throw std::invalid_argument("complex mode normalization dimensions mismatch");
    }

    std::complex<double> squaredNorm = 0.0;
    std::complex<double> slowness = halfSpaceSlowness(
        input.top, input, matrix, eigenvalue, rawMode[0]);

    int global = 0;
    for (std::size_t medium = firstAcoustic; medium <= lastAcoustic; ++medium)
    {
        const int count = matrix.counts[medium];
        const int offset = matrix.offsets[medium];
        const double spacing = matrix.spacing[medium];
        const double rhoMedium = matrix.rho[offset];
        const double rhoOmegaH2 =
            rhoMedium * matrix.omega2 * spacing * spacing;

        for (int local = 0; local <= count; ++local)
        {
            const double weight = (local == 0 || local == count) ? 0.5 : 1.0;
            const std::complex<double> phiSquared =
                rawMode[global + local] * rawMode[global + local];
            squaredNorm += weight * spacing * phiSquared / rhoMedium;
            slowness += weight * spacing *
                        (matrix.b1[offset + local] + 2.0) * phiSquared /
                        rhoOmegaH2;
        }
        global += count;
    }

    slowness += halfSpaceSlowness(
        input.bottom, input, matrix, eigenvalue,
        rawMode[rawMode.size() - 1]);
    const std::complex<double> topDerivative = boundaryAdmittanceDerivative(
        input.top, true, input, matrix, eigenvalue);
    const std::complex<double> bottomDerivative = boundaryAdmittanceDerivative(
        input.bottom, false, input, matrix, eigenvalue);
    const std::complex<double> renormalization =
        squaredNorm - topDerivative * rawMode[0] * rawMode[0] +
        bottomDerivative * rawMode[rawMode.size() - 1] *
            rawMode[rawMode.size() - 1];
    if (std::abs(renormalization) == 0.0)
    {
        throw std::runtime_error("complex mode has zero normalization integral");
    }

    std::complex<double> scale = 1.0 / std::sqrt(renormalization);
    if ((scale * rawMode[turningPoint]).real() < 0.0)
    {
        scale = -scale;
    }

    NormalizedComplexMode result;
    result.mode = scale * rawMode;
    result.complexNorm = renormalization * scale * scale;
    result.turningPoint = turningPoint;
    slowness *= scale * scale * matrix.omega / std::sqrt(eigenvalue);
    result.groupVelocity = (1.0 / slowness).real();
    result.scatterPerturbation = interfacialScatterPerturbation(
        input, matrix, eigenvalue, result.mode);
    if (!result.mode.allFinite() || !std::isfinite(result.groupVelocity))
    {
        throw std::runtime_error("complex mode normalization produced a non-finite result");
    }
    return result;
}
}
