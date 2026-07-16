#ifndef OPEN_OCEAN_KRAKENC_REFLECTION_BOUNDARY_H
#define OPEN_OCEAN_KRAKENC_REFLECTION_BOUNDARY_H

#include "algorithm/AcousticCase.h"
#include "algorithm/ComplexDispersion.h"

#include <complex>

namespace OpenOceanKrakenc
{
Impedance reflectionBoundaryImpedance(
    const AcousticBoundary &boundary,
    bool top,
    double omega2,
    std::complex<double> eigenvalue,
    std::complex<double> insideSoundSpeed);
}

#endif
