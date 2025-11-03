#ifndef SCATTER_H
#define SCATTER_H

#include "kkc_params.h"

std::complex<double> ScatterRoot(const std::complex<double>& z);
std::complex<double> KupIng(double sigma,
                           const std::complex<double>& eta1Sq,
                           double rho1,
                           const std::complex<double>& eta2Sq,
                           double rho2,
                           const std::complex<double>& P,
                           const std::complex<double>& U);

#endif // SCATTER_H