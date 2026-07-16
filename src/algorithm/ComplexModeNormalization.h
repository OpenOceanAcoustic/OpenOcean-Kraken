#ifndef OPEN_OCEAN_KRAKENC_COMPLEX_MODE_NORMALIZATION_H
#define OPEN_OCEAN_KRAKENC_COMPLEX_MODE_NORMALIZATION_H

#include "algorithm/AcousticCase.h"
#include "algorithm/ComplexMatrixBuilder.h"

#include <Eigen/Dense>

#include <complex>

namespace OpenOceanKrakenc
{
struct NormalizedComplexMode
{
    Eigen::VectorXcd mode;
    std::complex<double> complexNorm{};
    double groupVelocity = 0.0;
    int turningPoint = 0;
};

NormalizedComplexMode normalizeAcousticMode(
    const AcousticCase &input,
    const AcousticMatrix &matrix,
    std::complex<double> eigenvalue,
    int turningPoint,
    const Eigen::Ref<const Eigen::VectorXcd> &rawMode);
}

#endif
