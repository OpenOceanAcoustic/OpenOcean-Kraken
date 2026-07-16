#include "algorithm/InverseIterationComplex.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>

namespace OpenOceanKrakenc
{
namespace
{
bool finite(const Eigen::VectorXcd &values)
{
    for (Eigen::Index index = 0; index < values.size(); ++index)
    {
        if (!std::isfinite(values[index].real()) ||
            !std::isfinite(values[index].imag()))
        {
            return false;
        }
    }
    return true;
}

void multiplyTridiagonal(
    const Eigen::Ref<const Eigen::VectorXcd> &diagonal,
    const Eigen::Ref<const Eigen::VectorXcd> &offDiagonal,
    const Eigen::Ref<const Eigen::VectorXcd> &vector,
    Eigen::Ref<Eigen::VectorXcd> result)
{
    result = diagonal.array() * vector.array();
    for (Eigen::Index index = 0; index < offDiagonal.size(); ++index)
    {
        result[index] += offDiagonal[index] * vector[index + 1];
        result[index + 1] += offDiagonal[index] * vector[index];
    }
}
}

InverseIterationResult inverseIterationComplex(
    const Eigen::Ref<const Eigen::VectorXcd> &diagonal,
    const Eigen::Ref<const Eigen::VectorXcd> &offDiagonal,
    int maxIterations)
{
    InverseIterationResult result;
    if (diagonal.size() < 2 || offDiagonal.size() + 1 != diagonal.size())
    {
        result.failure = InverseIterationFailure::InvalidDimensions;
        return result;
    }
    if (maxIterations < 1)
    {
        result.failure = InverseIterationFailure::IterationLimit;
        return result;
    }
    if (!finite(diagonal) || !finite(offDiagonal))
    {
        result.failure = InverseIterationFailure::NonFiniteValue;
        return result;
    }

    const double matrixScale = std::max(
        1.0, diagonal.cwiseAbs().sum() + 2.0 * offDiagonal.cwiseAbs().sum());
    const double pivotFloor =
        100.0 * std::numeric_limits<double>::epsilon() * matrixScale;
    Eigen::VectorXcd iterate = Eigen::VectorXcd::Ones(diagonal.size());
    iterate.normalize();
    Eigen::VectorXcd pivots(diagonal.size());
    Eigen::VectorXcd rhs(diagonal.size());
    Eigen::VectorXcd next(diagonal.size());
    Eigen::VectorXcd residualVector(diagonal.size());

    for (int iteration = 1; iteration <= maxIterations; ++iteration)
    {
        pivots = diagonal;
        rhs = iterate;
        for (Eigen::Index row = 1; row < diagonal.size(); ++row)
        {
            if (std::abs(pivots[row - 1]) < pivotFloor)
            {
                pivots[row - 1] = std::complex<double>(pivotFloor, 0.0);
            }
            const std::complex<double> factor = offDiagonal[row - 1] / pivots[row - 1];
            pivots[row] -= factor * offDiagonal[row - 1];
            rhs[row] -= factor * rhs[row - 1];
        }
        if (std::abs(pivots[pivots.size() - 1]) < pivotFloor)
        {
            pivots[pivots.size() - 1] = std::complex<double>(pivotFloor, 0.0);
        }

        next[next.size() - 1] = rhs[rhs.size() - 1] / pivots[pivots.size() - 1];
        for (Eigen::Index row = diagonal.size() - 1; row-- > 0;)
        {
            next[row] = (rhs[row] - offDiagonal[row] * next[row + 1]) / pivots[row];
        }
        const double nextNorm = next.norm();
        if (!finite(next) || nextNorm == 0.0)
        {
            result.failure = InverseIterationFailure::NonFiniteValue;
            result.iterations = iteration;
            return result;
        }
        next /= nextNorm;

        Eigen::Index phaseIndex = 0;
        next.cwiseAbs().maxCoeff(&phaseIndex);
        const std::complex<double> anchor = next[phaseIndex];
        next *= std::conj(anchor) / std::abs(anchor);

        multiplyTridiagonal(diagonal, offDiagonal, next, residualVector);
        const double residual = residualVector.norm();
        result.vector = next;
        result.iterations = iteration;
        result.residual = residual;
        if (residual <= 1.0e-10 * matrixScale)
        {
            result.converged = true;
            result.failure = InverseIterationFailure::None;
            return result;
        }
        iterate = next;
    }

    result.failure = InverseIterationFailure::IterationLimit;
    return result;
}
}
