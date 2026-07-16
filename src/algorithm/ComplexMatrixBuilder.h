#ifndef OPEN_OCEAN_KRAKENC_COMPLEX_MATRIX_BUILDER_H
#define OPEN_OCEAN_KRAKENC_COMPLEX_MATRIX_BUILDER_H

#include "algorithm/AcousticCase.h"

#include <complex>
#include <vector>

namespace OpenOceanKrakenc
{
struct AcousticMatrix
{
    double omega = 0.0;
    double omega2 = 0.0;
    std::vector<int> counts;
    std::vector<int> offsets;
    std::vector<double> spacing;
    std::vector<bool> layerElastic;
    std::vector<double> depth;
    std::vector<double> rho;
    std::vector<std::complex<double>> cp;
    std::vector<std::complex<double>> cs;
    std::vector<std::complex<double>> b1;
    std::vector<std::complex<double>> b2;
    std::vector<std::complex<double>> b3;
    std::vector<std::complex<double>> b4;
    std::vector<std::complex<double>> dynamicRho;
};

AcousticMatrix buildAcousticMatrix(const AcousticCase &input, int meshMultiplier);
}

#endif
