#ifndef OPEN_OCEAN_KRAKENC_SSP_INTERPOLATION_H
#define OPEN_OCEAN_KRAKENC_SSP_INTERPOLATION_H

#include "algorithm/AcousticCase.h"

#include <array>
#include <complex>
#include <cstddef>
#include <vector>

namespace OpenOceanKrakenc
{
struct InterpolatedMaterial
{
    std::complex<double> cp{};
    std::complex<double> cs{};
    double rho = 0.0;
};

void validateSspLayer(const AcousticLayer &layer,
                      std::size_t mediumNumber);

class PreparedComplexSspLayer
{
public:
    PreparedComplexSspLayer(const AcousticLayer &layer,
                            AcousticInterpolation interpolation,
                            const AttenuationContext &context,
                            std::size_t mediumNumber);

    InterpolatedMaterial evaluate(double depth) const;
    const InterpolatedMaterial &top() const noexcept;
    const InterpolatedMaterial &bottom() const noexcept;
    bool elastic() const noexcept;

private:
    using SegmentCoefficients =
        std::array<std::complex<double>, 4>;

    AcousticInterpolation interpolation_;
    std::size_t mediumNumber_;
    bool elastic_ = false;
    std::vector<double> depths_;
    std::vector<InterpolatedMaterial> nodes_;
    std::vector<SegmentCoefficients> cpCoefficients_;
    std::vector<SegmentCoefficients> csCoefficients_;
    std::vector<SegmentCoefficients> rhoCoefficients_;
};

double interpolateRealCp(const AcousticLayer &layer,
                         AcousticInterpolation interpolation,
                         double depth,
                         std::size_t mediumNumber);
}

#endif
