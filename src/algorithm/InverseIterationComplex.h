#ifndef OPEN_OCEAN_KRAKENC_INVERSE_ITERATION_COMPLEX_H
#define OPEN_OCEAN_KRAKENC_INVERSE_ITERATION_COMPLEX_H

#include <Eigen/Dense>

namespace OpenOceanKrakenc
{
enum class InverseIterationFailure
{
    None,
    InvalidDimensions,
    NonFiniteValue,
    IterationLimit
};

struct InverseIterationResult
{
    Eigen::VectorXcd vector;
    int iterations = 0;
    double residual = 0.0;
    bool converged = false;
    InverseIterationFailure failure = InverseIterationFailure::None;
};

InverseIterationResult inverseIterationComplex(
    const Eigen::Ref<const Eigen::VectorXcd> &diagonal,
    const Eigen::Ref<const Eigen::VectorXcd> &offDiagonal,
    int maxIterations = 30);
}

#endif
