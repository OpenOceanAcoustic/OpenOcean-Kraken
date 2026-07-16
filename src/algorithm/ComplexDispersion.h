#ifndef OPEN_OCEAN_KRAKENC_COMPLEX_DISPERSION_H
#define OPEN_OCEAN_KRAKENC_COMPLEX_DISPERSION_H

#include "algorithm/AcousticCase.h"
#include "algorithm/ComplexMatrixBuilder.h"
#include "algorithm/ComplexRootFinder.h"

#include <complex>
#include <vector>

namespace OpenOceanKrakenc
{
struct Impedance
{
    std::complex<double> f{};
    std::complex<double> g{};
    int power10 = 0;
};

Impedance acousticBoundaryImpedance(const AcousticBoundary &boundary,
                                    bool top,
                                    double omega2,
                                    std::complex<double> eigenvalue,
                                    char attenuationUnit,
                                    double frequency);

class ComplexDispersion
{
public:
    ComplexDispersion(const AcousticCase &input, const AcousticMatrix &matrix);

    ScaledComplex evaluate(
        std::complex<double> eigenvalue,
        const std::vector<std::complex<double>> &acceptedRoots) const;

private:
    Impedance boundaryImpedance(bool top,
                                std::complex<double> eigenvalue) const;

    const AcousticCase &input_;
    const AcousticMatrix &matrix_;
    std::size_t firstAcoustic_ = 0;
    std::size_t lastAcoustic_ = 0;
    std::complex<double> topCp_{};
    std::complex<double> topCs_{};
    std::complex<double> bottomCp_{};
    std::complex<double> bottomCs_{};
};
}

#endif
