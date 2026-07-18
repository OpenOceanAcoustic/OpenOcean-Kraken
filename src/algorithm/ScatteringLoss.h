#ifndef OPEN_OCEAN_KRAKENC_SCATTERING_LOSS_H
#define OPEN_OCEAN_KRAKENC_SCATTERING_LOSS_H

#include <complex>

namespace OpenOceanKrakenc
{
std::complex<double> scatterRoot(std::complex<double> value);

std::complex<double> kupermanIngenito(
    double sigma,
    std::complex<double> eta1Squared,
    double rho1,
    std::complex<double> eta2Squared,
    double rho2,
    std::complex<double> pressure,
    std::complex<double> pressureDerivativeOverDensity);
}

#endif
