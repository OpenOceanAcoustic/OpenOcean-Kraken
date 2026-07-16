#ifndef OPEN_OCEAN_KRAKENC_COMPLEX_NUMERICS_H
#define OPEN_OCEAN_KRAKENC_COMPLEX_NUMERICS_H

#include "OpenOceanKrakencParams.h"

#include <complex>

namespace OpenOceanKrakenc
{
std::complex<double> pekerisRoot(std::complex<double> value);

std::complex<double> complexSoundSpeed(double soundSpeed,
                                       double attenuation,
                                       double frequency,
                                       char attenuationUnit);
}

#endif
