#ifndef OPEN_OCEAN_KRAKENC_ELASTIC_COMPOUND_H
#define OPEN_OCEAN_KRAKENC_ELASTIC_COMPOUND_H

#include "algorithm/AcousticCase.h"
#include "algorithm/ComplexDispersion.h"
#include "algorithm/ComplexMatrixBuilder.h"

#include <array>
#include <complex>

namespace OpenOceanKrakenc
{
using CompoundState = std::array<std::complex<double>, 5>;

bool compoundStateFinite(const CompoundState &state);

void propagateElasticDown(const AcousticMatrix &matrix,
                          std::size_t layer,
                          std::complex<double> eigenvalue,
                          CompoundState &state,
                          int &power10);

void propagateElasticUp(const AcousticMatrix &matrix,
                        std::size_t layer,
                        std::complex<double> eigenvalue,
                        CompoundState &state,
                        int &power10);

Impedance caseBoundaryImpedance(const AcousticCase &input,
                                const AcousticMatrix &matrix,
                                bool top,
                                std::complex<double> eigenvalue);
}

#endif
