#ifndef OPEN_OCEAN_KRAKENC_COMPLEX_MODE_SOLVER_H
#define OPEN_OCEAN_KRAKENC_COMPLEX_MODE_SOLVER_H

#include "algorithm/AcousticCase.h"
#include "algorithm/ComplexMatrixBuilder.h"
#include "algorithm/InverseIterationComplex.h"

#include <Eigen/Dense>

#include <complex>

namespace OpenOceanKrakenc
{
struct AcousticModeMatrix
{
    Eigen::VectorXd depth;
    Eigen::VectorXcd diagonal;
    Eigen::VectorXcd offDiagonal;
    int turningPoint = 0;
};

struct ComplexModeResult
{
    Eigen::VectorXd depth;
    Eigen::VectorXcd mode;
    int turningPoint = 0;
    int iterations = 0;
    double residual = 0.0;
    bool converged = false;
    InverseIterationFailure failure = InverseIterationFailure::None;
};

AcousticModeMatrix buildAcousticModeMatrix(
    const AcousticCase &input,
    const AcousticMatrix &matrix,
    std::complex<double> eigenvalue);

ComplexModeResult solveAcousticMode(
    const AcousticCase &input,
    const AcousticMatrix &matrix,
    std::complex<double> eigenvalue,
    int maxIterations = 30);
}

#endif
