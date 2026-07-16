#include "algorithm/ComplexModeSolver.h"

#include "algorithm/ComplexDispersion.h"
#include "algorithm/ElasticCompound.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace OpenOceanKrakenc
{
namespace
{
void validateMetadata(const AcousticCase &input, const AcousticMatrix &matrix)
{
    const std::size_t media = input.layers.size();
    if (media == 0 || matrix.counts.size() != media ||
        matrix.offsets.size() != media || matrix.spacing.size() != media)
    {
        throw std::invalid_argument("acoustic mode matrix metadata mismatch");
    }
}
}

AcousticModeMatrix buildAcousticModeMatrix(
    const AcousticCase &input,
    const AcousticMatrix &matrix,
    std::complex<double> eigenvalue)
{
    validateMetadata(input, matrix);

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
    if (firstAcoustic == input.layers.size())
    {
        throw std::invalid_argument("acoustic mode extraction requires an acoustic layer");
    }
    int intervalCount = 0;
    for (std::size_t medium = firstAcoustic; medium <= lastAcoustic; ++medium)
    {
        if (matrix.layerElastic[medium])
        {
            throw std::invalid_argument("acoustic layers must form one contiguous block");
        }
        intervalCount += matrix.counts[medium];
    }

    AcousticModeMatrix result;
    result.depth.resize(intervalCount + 1);
    result.diagonal = Eigen::VectorXcd::Zero(intervalCount + 1);
    result.offDiagonal = Eigen::VectorXcd::Zero(intervalCount);
    result.turningPoint = intervalCount;

    int global = 0;
    for (std::size_t medium = firstAcoustic; medium <= lastAcoustic; ++medium)
    {
        const int count = matrix.counts[medium];
        const int offset = matrix.offsets[medium];
        const double spacing = matrix.spacing[medium];
        const double hRho = spacing * matrix.rho[offset];
        if (!(hRho > 0.0))
        {
            throw std::invalid_argument("acoustic mode matrix requires positive density");
        }
        const std::complex<double> xh2 = eigenvalue * spacing * spacing;

        for (int local = 0; local <= count; ++local)
        {
            const int node = global + local;
            const std::complex<double> value =
                (matrix.b1[offset + local] - xh2) / hRho;
            result.depth[node] = matrix.depth[offset + local];
            if (medium > firstAcoustic && local == 0)
            {
                result.diagonal[node] = 0.5 * (result.diagonal[node] + value);
            }
            else
            {
                result.diagonal[node] = value;
            }
            if (local < count)
            {
                result.offDiagonal[node] = 1.0 / hRho;
            }
            if (local > 0 &&
                (matrix.b1[offset + local] - xh2).real() + 2.0 > 0.0)
            {
                result.turningPoint = std::min(result.turningPoint, node);
            }
        }
        global += count;
    }

    const Impedance top = caseBoundaryImpedance(
        input, matrix, true, eigenvalue);
    if (top.g == std::complex<double>(0.0, 0.0))
    {
        result.diagonal[0] = 1.0;
        result.offDiagonal[0] = 0.0;
    }
    else
    {
        result.diagonal[0] = 0.5 * result.diagonal[0] + top.f / top.g;
    }

    const Impedance bottom = caseBoundaryImpedance(
        input, matrix, false, eigenvalue);
    const Eigen::Index last = result.diagonal.size() - 1;
    if (bottom.g == std::complex<double>(0.0, 0.0))
    {
        result.diagonal[last] = 1.0;
        result.offDiagonal[last - 1] = 0.0;
    }
    else
    {
        result.diagonal[last] =
            0.5 * result.diagonal[last] - bottom.f / bottom.g;
    }
    return result;
}

ComplexModeResult solveAcousticMode(
    const AcousticCase &input,
    const AcousticMatrix &matrix,
    std::complex<double> eigenvalue,
    int maxIterations)
{
    const AcousticModeMatrix assembled =
        buildAcousticModeMatrix(input, matrix, eigenvalue);
    const InverseIterationResult inverse = inverseIterationComplex(
        assembled.diagonal, assembled.offDiagonal, maxIterations);

    ComplexModeResult result;
    result.depth = assembled.depth;
    result.mode = inverse.vector;
    result.turningPoint = assembled.turningPoint;
    result.iterations = inverse.iterations;
    result.residual = inverse.residual;
    result.converged = inverse.converged;
    result.failure = inverse.failure;
    return result;
}
}
