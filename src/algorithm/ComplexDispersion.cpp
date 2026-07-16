#include "algorithm/ComplexDispersion.h"

#include "algorithm/ComplexNumerics.h"
#include "algorithm/ElasticCompound.h"
#include "algorithm/ReflectionBoundary.h"

#include <cmath>
#include <stdexcept>

namespace OpenOceanKrakenc
{
namespace
{
constexpr double roof = 1.0e50;
constexpr double floorValue = 1.0e-50;
constexpr int scalePower = 50;

void scaleRecurrence(std::complex<double> &p0,
                     std::complex<double> &p1,
                     std::complex<double> &p2,
                     int &power10)
{
    while (std::abs(p2.real()) > roof)
    {
        p0 *= floorValue;
        p1 *= floorValue;
        p2 *= floorValue;
        power10 += scalePower;
    }
}

void scaleValue(std::complex<double> &value, int &power10)
{
    while (std::abs(value.real()) > roof)
    {
        value *= floorValue;
        power10 += scalePower;
    }
    while (std::abs(value.real()) < floorValue &&
           value != std::complex<double>(0.0, 0.0))
    {
        value *= roof;
        power10 -= scalePower;
    }
}
}

Impedance acousticBoundaryImpedance(const AcousticBoundary &boundary,
                                    bool top,
                                    double omega2,
                                    std::complex<double> eigenvalue,
                                    char attenuationUnit,
                                    double frequency)
{
    Impedance result;
    switch (boundary.type)
    {
    case AcousticBoundaryType::Vacuum:
        result.f = 1.0;
        result.g = 0.0;
        break;
    case AcousticBoundaryType::Rigid:
        result.f = 0.0;
        result.g = 1.0;
        break;
    case AcousticBoundaryType::HalfSpace:
    {
        if (!(boundary.cp > 0.0) || !(boundary.rho > 0.0))
        {
            throw std::invalid_argument("half-space must have positive cp and rho");
        }
        const std::complex<double> cp = complexSoundSpeed(
            boundary.cp, boundary.alphaP, frequency, attenuationUnit);
        if (boundary.cs <= 0.0)
        {
            result.f = pekerisRoot(eigenvalue - omega2 / (cp * cp));
            result.g = boundary.rho;
        }
        else
        {
            const std::complex<double> cs = complexSoundSpeed(
                boundary.cs, boundary.alphaS, frequency, attenuationUnit);
            const std::complex<double> gammaS2 = eigenvalue - omega2 / (cs * cs);
            const std::complex<double> gammaP2 = eigenvalue - omega2 / (cp * cp);
            const std::complex<double> gammaS = pekerisRoot(gammaS2);
            const std::complex<double> gammaP = pekerisRoot(gammaP2);
            const std::complex<double> mu = boundary.rho * cs * cs;
            result.f = omega2 * gammaP * (eigenvalue - gammaS2);
            result.g = ((gammaS2 + eigenvalue) * (gammaS2 + eigenvalue) -
                        4.0 * gammaS * gammaP * eigenvalue) * mu;
        }
        break;
    }
    case AcousticBoundaryType::ReflectionCoefficient:
    case AcousticBoundaryType::InternalReflection:
        throw std::invalid_argument("reflection-table impedance requires the interior sound speed");
    }
    if (top)
    {
        result.g = -result.g;
    }
    return result;
}

ComplexDispersion::ComplexDispersion(const AcousticCase &input,
                                     const AcousticMatrix &matrix)
    : input_(input), matrix_(matrix)
{
    if (matrix_.counts.size() != input_.layers.size() ||
        matrix_.offsets.size() != input_.layers.size() ||
        matrix_.spacing.size() != input_.layers.size() ||
        matrix_.layerElastic.size() != input_.layers.size())
    {
        throw std::invalid_argument("acoustic matrix metadata does not match the case");
    }
    firstAcoustic_ = input_.layers.size();
    lastAcoustic_ = input_.layers.size();
    for (std::size_t medium = 0; medium < input_.layers.size(); ++medium)
    {
        if (!matrix_.layerElastic[medium])
        {
            if (firstAcoustic_ == input_.layers.size())
            {
                firstAcoustic_ = medium;
            }
            lastAcoustic_ = medium;
        }
    }
    if (firstAcoustic_ == input_.layers.size())
    {
        throw std::invalid_argument("complex dispersion requires an acoustic layer");
    }
    for (std::size_t medium = firstAcoustic_; medium <= lastAcoustic_; ++medium)
    {
        if (matrix_.layerElastic[medium])
        {
            throw std::invalid_argument("acoustic layers must form one contiguous block");
        }
    }
    const auto cacheSpeeds = [&](const AcousticBoundary &boundary,
                                 std::complex<double> &cp,
                                 std::complex<double> &cs) {
        if (boundary.cp > 0.0)
        {
            cp = complexSoundSpeed(boundary.cp, boundary.alphaP,
                                   input_.frequency, input_.attenuationUnit);
        }
        if (boundary.cs > 0.0)
        {
            cs = complexSoundSpeed(boundary.cs, boundary.alphaS,
                                   input_.frequency, input_.attenuationUnit);
        }
    };
    cacheSpeeds(input_.top, topCp_, topCs_);
    cacheSpeeds(input_.bottom, bottomCp_, bottomCs_);
}

Impedance ComplexDispersion::boundaryImpedance(
    bool top,
    std::complex<double> eigenvalue) const
{
    const bool direct = top ? firstAcoustic_ == 0
                            : lastAcoustic_ + 1 == matrix_.layerElastic.size();
    if (!direct)
    {
        return caseBoundaryImpedance(input_, matrix_, top, eigenvalue);
    }
    const AcousticBoundary &boundary = top ? input_.top : input_.bottom;
    if (boundary.type == AcousticBoundaryType::ReflectionCoefficient ||
        boundary.type == AcousticBoundaryType::InternalReflection)
    {
        const int interiorIndex = top
            ? matrix_.offsets[firstAcoustic_]
            : matrix_.offsets[lastAcoustic_] + matrix_.counts[lastAcoustic_];
        return reflectionBoundaryImpedance(
            boundary, top, matrix_.omega2, eigenvalue,
            matrix_.cp[interiorIndex]);
    }
    Impedance result;
    if (boundary.type == AcousticBoundaryType::Vacuum)
    {
        result.f = 1.0;
        return result;
    }
    if (boundary.type == AcousticBoundaryType::Rigid)
    {
        result.g = top ? -1.0 : 1.0;
        return result;
    }
    const std::complex<double> cp = top ? topCp_ : bottomCp_;
    const std::complex<double> cs = top ? topCs_ : bottomCs_;
    if (cs == std::complex<double>(0.0, 0.0))
    {
        result.f = pekerisRoot(eigenvalue - matrix_.omega2 / (cp * cp));
        result.g = boundary.rho;
    }
    else
    {
        const std::complex<double> gammaS2 =
            eigenvalue - matrix_.omega2 / (cs * cs);
        const std::complex<double> gammaP2 =
            eigenvalue - matrix_.omega2 / (cp * cp);
        const std::complex<double> gammaS = pekerisRoot(gammaS2);
        const std::complex<double> gammaP = pekerisRoot(gammaP2);
        const std::complex<double> mu = boundary.rho * cs * cs;
        result.f = matrix_.omega2 * gammaP * (eigenvalue - gammaS2);
        result.g = ((gammaS2 + eigenvalue) * (gammaS2 + eigenvalue) -
                    4.0 * gammaS * gammaP * eigenvalue) * mu;
    }
    if (top)
    {
        result.g = -result.g;
    }
    return result;
}

ScaledComplex ComplexDispersion::evaluate(
    std::complex<double> eigenvalue,
    const std::vector<std::complex<double>> &acceptedRoots) const
{
    Impedance running = boundaryImpedance(false, eigenvalue);

    for (std::size_t reverse = lastAcoustic_ + 1; reverse-- > firstAcoustic_;)
    {
        const int count = matrix_.counts[reverse];
        const int offset = matrix_.offsets[reverse];
        const double spacing = matrix_.spacing[reverse];
        const std::complex<double> h2k2 = spacing * spacing * eigenvalue;
        const double rhoMedium = matrix_.rho[offset];

        std::complex<double> p0 = 0.0;
        std::complex<double> p1 = -2.0 * running.g;
        std::complex<double> p2 =
            (matrix_.b1[offset + count] - h2k2) * running.g -
            2.0 * spacing * running.f * rhoMedium;

        for (int index = offset + count - 1; index >= offset; --index)
        {
            p0 = p1;
            p1 = p2;
            p2 = (h2k2 - matrix_.b1[index]) * p1 - p0;
            scaleRecurrence(p0, p1, p2, running.power10);
        }

        running.f = -(p2 - p0) / (2.0 * spacing * rhoMedium);
        running.g = -p1;
    }

    const Impedance top = boundaryImpedance(true, eigenvalue);
    std::complex<double> determinant = running.f * top.g - running.g * top.f;
    int power10 = running.power10 + top.power10;
    scaleValue(determinant, power10);

    std::size_t rootsSinceScaling = 0;
    for (const std::complex<double> root : acceptedRoots)
    {
        const std::complex<double> difference = eigenvalue - root;
        if (difference == std::complex<double>(0.0, 0.0))
        {
            return {{std::numeric_limits<double>::infinity(), 0.0}, 0};
        }
        determinant /= difference;
        if (++rootsSinceScaling == 8)
        {
            scaleValue(determinant, power10);
            rootsSinceScaling = 0;
        }
    }
    scaleValue(determinant, power10);
    return {determinant, power10};
}
}
