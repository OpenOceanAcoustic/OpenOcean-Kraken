#include "algorithm/ComplexModeNormalization.h"

#include "algorithm/ComplexDispersion.h"
#include "algorithm/ComplexNumerics.h"
#include "algorithm/ElasticCompound.h"

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
        boundary.cp, boundary.alphaP, input.frequency,
        input.attenuationUnit);
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
    if (!result.mode.allFinite() || !std::isfinite(result.groupVelocity))
    {
        throw std::runtime_error("complex mode normalization produced a non-finite result");
    }
    return result;
}
}
