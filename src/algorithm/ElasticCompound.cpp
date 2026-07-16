#include "algorithm/ElasticCompound.h"

#include "algorithm/ComplexNumerics.h"
#include "algorithm/ReflectionBoundary.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace OpenOceanKrakenc
{
namespace
{
constexpr double roof = 1.0e50;
constexpr double floorValue = 1.0e-50;

CompoundState derivative(const AcousticMatrix &matrix,
                         int index,
                         std::complex<double> eigenvalue,
                         const CompoundState &state)
{
    const std::complex<double> twoX = 2.0 * eigenvalue;
    const double twoH = 2.0 * matrix.spacing[
        static_cast<std::size_t>(std::upper_bound(
            matrix.offsets.begin(), matrix.offsets.end(), index) -
            matrix.offsets.begin() - 1)];
    const std::complex<double> fourHX = 2.0 * twoH * eigenvalue;
    const std::complex<double> xB3 =
        eigenvalue * matrix.b3[index] - matrix.dynamicRho[index];
    return {{
        matrix.b1[index] * state[3] - matrix.b2[index] * state[4],
        -matrix.dynamicRho[index] * state[3] - xB3 * state[4],
        twoH * state[3] + matrix.b4[index] * state[4],
        xB3 * state[0] + matrix.b2[index] * state[1] -
            twoX * matrix.b4[index] * state[2],
        matrix.dynamicRho[index] * state[0] - matrix.b1[index] * state[1] -
            fourHX * state[2],
    }};
}

CompoundState addScaled(const CompoundState &base,
                        const CompoundState &increment,
                        double scale)
{
    CompoundState result;
    for (std::size_t index = 0; index < result.size(); ++index)
    {
        result[index] = base[index] + scale * increment[index];
    }
    return result;
}

void rescaleIfNeeded(CompoundState &current, CompoundState &previous, int &power10)
{
    const double pivot = std::abs(current[1].real());
    double scale = 1.0;
    if (pivot < floorValue)
    {
        scale = roof;
        power10 -= 50;
    }
    else if (pivot > roof)
    {
        scale = floorValue;
        power10 += 50;
    }
    if (scale != 1.0)
    {
        for (std::size_t index = 0; index < current.size(); ++index)
        {
            current[index] *= scale;
            previous[index] *= scale;
        }
    }
}

CompoundState boundaryState(const AcousticBoundary &boundary,
                            const AcousticCase &input,
                            const AcousticMatrix &matrix,
                            std::complex<double> eigenvalue)
{
    CompoundState state{};
    if (boundary.type == AcousticBoundaryType::Vacuum)
    {
        state[0] = 1.0;
        return state;
    }
    if (boundary.type == AcousticBoundaryType::Rigid)
    {
        state[1] = 1.0;
        return state;
    }

    const std::complex<double> cp = complexSoundSpeed(
        boundary.cp, boundary.alphaP, input.frequency, input.attenuationUnit);
    if (boundary.cs <= 0.0)
    {
        state[0] = pekerisRoot(eigenvalue - matrix.omega2 / (cp * cp));
        state[1] = boundary.rho;
        return state;
    }
    const std::complex<double> cs = complexSoundSpeed(
        boundary.cs, boundary.alphaS, input.frequency, input.attenuationUnit);
    const std::complex<double> gammaS2 =
        eigenvalue - matrix.omega2 / (cs * cs);
    const std::complex<double> gammaP2 =
        eigenvalue - matrix.omega2 / (cp * cp);
    const std::complex<double> gammaS = pekerisRoot(gammaS2);
    const std::complex<double> gammaP = pekerisRoot(gammaP2);
    const std::complex<double> mu = boundary.rho * cs * cs;
    state[0] = (gammaS * gammaP - eigenvalue) / mu;
    state[1] = ((gammaS2 + eigenvalue) * (gammaS2 + eigenvalue) -
                4.0 * gammaS * gammaP * eigenvalue) * mu;
    state[2] = 2.0 * gammaS * gammaP - gammaS2 - eigenvalue;
    state[3] = gammaP * (eigenvalue - gammaS2);
    state[4] = gammaS * (gammaS2 - eigenvalue);
    return state;
}
}

bool compoundStateFinite(const CompoundState &state)
{
    for (const std::complex<double> value : state)
    {
        if (!std::isfinite(value.real()) || !std::isfinite(value.imag()))
        {
            return false;
        }
    }
    return true;
}

void propagateElasticDown(const AcousticMatrix &matrix,
                          std::size_t layer,
                          std::complex<double> eigenvalue,
                          CompoundState &state,
                          int &power10)
{
    if (layer >= matrix.layerElastic.size() || !matrix.layerElastic[layer])
    {
        throw std::invalid_argument("downward compound propagation requires an elastic layer");
    }
    const int offset = matrix.offsets[layer];
    const int count = matrix.counts[layer];
    CompoundState previous = state;
    CompoundState current = addScaled(state, derivative(matrix, offset, eigenvalue, state), 0.5);
    CompoundState next{};
    for (int local = 1; local <= count; ++local)
    {
        previous = state;
        state = current;
        current = addScaled(previous,
                            derivative(matrix, offset + local, eigenvalue, state), 1.0);
        if (local != count)
        {
            rescaleIfNeeded(current, state, power10);
        }
    }
    for (std::size_t index = 0; index < state.size(); ++index)
    {
        next[index] = (previous[index] + 2.0 * state[index] + current[index]) / 4.0;
    }
    state = next;
}

void propagateElasticUp(const AcousticMatrix &matrix,
                        std::size_t layer,
                        std::complex<double> eigenvalue,
                        CompoundState &state,
                        int &power10)
{
    if (layer >= matrix.layerElastic.size() || !matrix.layerElastic[layer])
    {
        throw std::invalid_argument("upward compound propagation requires an elastic layer");
    }
    const int offset = matrix.offsets[layer];
    const int count = matrix.counts[layer];
    CompoundState previous = state;
    CompoundState current = addScaled(
        state, derivative(matrix, offset + count, eigenvalue, state), -0.5);
    CompoundState next{};
    for (int local = count - 1; local >= 0; --local)
    {
        previous = state;
        state = current;
        current = addScaled(previous,
                            derivative(matrix, offset + local, eigenvalue, state), -1.0);
        if (local != 0)
        {
            rescaleIfNeeded(current, state, power10);
        }
    }
    for (std::size_t index = 0; index < state.size(); ++index)
    {
        next[index] = (previous[index] + 2.0 * state[index] + current[index]) / 4.0;
    }
    state = next;
}

Impedance caseBoundaryImpedance(const AcousticCase &input,
                                const AcousticMatrix &matrix,
                                bool top,
                                std::complex<double> eigenvalue)
{
    std::size_t firstAcoustic = matrix.layerElastic.size();
    std::size_t lastAcoustic = matrix.layerElastic.size();
    for (std::size_t layer = 0; layer < matrix.layerElastic.size(); ++layer)
    {
        if (!matrix.layerElastic[layer])
        {
            if (firstAcoustic == matrix.layerElastic.size())
            {
                firstAcoustic = layer;
            }
            lastAcoustic = layer;
        }
    }
    if (firstAcoustic == matrix.layerElastic.size())
    {
        throw std::invalid_argument("at least one acoustic layer is required");
    }
    if ((top && firstAcoustic == 0) ||
        (!top && lastAcoustic + 1 == matrix.layerElastic.size()))
    {
        const AcousticBoundary &boundary = top ? input.top : input.bottom;
        if (boundary.type == AcousticBoundaryType::ReflectionCoefficient ||
            boundary.type == AcousticBoundaryType::InternalReflection)
        {
            const int interiorIndex = top
                                          ? matrix.offsets[firstAcoustic]
                                          : matrix.offsets[lastAcoustic] +
                                                matrix.counts[lastAcoustic];
            return reflectionBoundaryImpedance(
                boundary, top, matrix.omega2, eigenvalue,
                matrix.cp[interiorIndex]);
        }
        return acousticBoundaryImpedance(
            boundary, top, matrix.omega2, eigenvalue,
            input.attenuationUnit, input.frequency);
    }

    CompoundState state = boundaryState(
        top ? input.top : input.bottom, input, matrix, eigenvalue);
    int power10 = 0;
    if (top)
    {
        for (std::size_t layer = 0; layer < firstAcoustic; ++layer)
        {
            propagateElasticDown(matrix, layer, eigenvalue, state, power10);
        }
    }
    else
    {
        for (std::size_t layer = matrix.layerElastic.size(); layer-- > lastAcoustic + 1;)
        {
            propagateElasticUp(matrix, layer, eigenvalue, state, power10);
        }
    }
    return {matrix.omega2 * state[3], state[1], power10};
}
}
