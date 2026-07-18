#include "algorithm/ScatteringLoss.h"

#include <cmath>
#include <stdexcept>

namespace OpenOceanKrakenc
{
std::complex<double> scatterRoot(std::complex<double> value)
{
    if (value.real() >= 0.0)
    {
        return std::sqrt(value);
    }
    return -std::complex<double>(0.0, 1.0) * std::sqrt(-value);
}

std::complex<double> kupermanIngenito(
    double sigma,
    std::complex<double> eta1Squared,
    double rho1,
    std::complex<double> eta2Squared,
    double rho2,
    std::complex<double> pressure,
    std::complex<double> pressureDerivativeOverDensity)
{
    if (sigma < 0.0 || !(rho1 > 0.0) || !(rho2 > 0.0))
    {
        throw std::invalid_argument(
            "Kuperman-Ingenito requires non-negative roughness and positive densities");
    }
    if (sigma == 0.0)
    {
        return 0.0;
    }
    const std::complex<double> imaginary(0.0, 1.0);
    const std::complex<double> eta1 = scatterRoot(eta1Squared);
    const std::complex<double> eta2 = scatterRoot(eta2Squared);
    const std::complex<double> denominator = rho1 * eta2 + rho2 * eta1;
    if (denominator == std::complex<double>(0.0, 0.0))
    {
        return 0.0;
    }
    const std::complex<double> a11 =
        0.5 * (eta1Squared - eta2Squared) -
        (rho2 * eta1Squared - rho1 * eta2Squared) * (eta1 + eta2) /
            denominator;
    const std::complex<double> a12 =
        imaginary * std::pow(rho2 - rho1, 2.0) * eta1 * eta2 /
        denominator;
    const std::complex<double> a21 =
        -imaginary * std::pow(rho2 * eta1Squared - rho1 * eta2Squared, 2.0) /
        (rho1 * rho2 * denominator);
    const std::complex<double> a22 =
        0.5 * (eta1Squared - eta2Squared) +
        (rho2 - rho1) * eta1 * eta2 * (eta1 + eta2) / denominator;
    return -sigma * sigma *
           (-a21 * pressure * pressure +
            (a11 - a22) * pressure * pressureDerivativeOverDensity +
            a12 * pressureDerivativeOverDensity *
                pressureDerivativeOverDensity);
}
}
