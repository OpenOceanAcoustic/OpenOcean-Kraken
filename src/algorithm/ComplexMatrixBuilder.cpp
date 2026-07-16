#include "algorithm/ComplexMatrixBuilder.h"

#include "algorithm/ComplexNumerics.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace OpenOceanKrakenc
{
AcousticMatrix buildAcousticMatrix(const AcousticCase &input, int meshMultiplier)
{
    if (meshMultiplier < 1)
    {
        throw std::invalid_argument("mesh multiplier must be positive");
    }
    if (!(input.frequency > 0.0) || input.layers.empty())
    {
        throw std::invalid_argument("acoustic case is incomplete");
    }

    AcousticMatrix matrix;
    matrix.omega = 2.0 * pi * input.frequency;
    matrix.omega2 = matrix.omega * matrix.omega;
    int offset = 0;

    for (const AcousticLayer &layer : input.layers)
    {
        const bool elastic = layer.samples.front().cs > 0.0;
        for (const AcousticSample &sample : layer.samples)
        {
            if ((sample.cs > 0.0) != elastic)
            {
                throw std::invalid_argument("a medium cannot mix acoustic and elastic samples");
            }
        }

        const int count = static_cast<int>(
            layer.baseMesh * meshMultiplier * input.frequency / input.referenceFrequency);
        if (count < 1)
        {
            throw std::invalid_argument("scaled acoustic mesh must contain at least one interval");
        }
        const double spacing = (layer.bottomDepth - layer.topDepth) / count;
        matrix.counts.push_back(count);
        matrix.offsets.push_back(offset);
        matrix.spacing.push_back(spacing);
        matrix.layerElastic.push_back(elastic);

        std::vector<std::complex<double>> sampleCp(layer.samples.size());
        std::vector<std::complex<double>> sampleCs(layer.samples.size(), 0.0);
        for (std::size_t sample = 0; sample < layer.samples.size(); ++sample)
        {
            sampleCp[sample] = complexSoundSpeed(
                layer.samples[sample].cp, layer.samples[sample].alphaP,
                input.frequency, input.attenuationUnit);
            if (elastic)
            {
                sampleCs[sample] = complexSoundSpeed(
                    layer.samples[sample].cs, layer.samples[sample].alphaS,
                    input.frequency, input.attenuationUnit);
            }
        }
        std::size_t lowerIndex = 0;

        for (int local = 0; local <= count; ++local)
        {
            double depth = layer.topDepth + local * spacing;
            if (local == count)
            {
                depth = layer.bottomDepth;
            }
            while (lowerIndex + 1 < layer.samples.size() - 1 &&
                   depth >= layer.samples[lowerIndex + 1].depth)
            {
                ++lowerIndex;
            }
            const AcousticSample &lower = layer.samples[lowerIndex];
            const AcousticSample &upper = layer.samples[lowerIndex + 1];
            const double fraction = (depth - lower.depth) / (upper.depth - lower.depth);

            const std::complex<double> lowerCp = sampleCp[lowerIndex];
            const std::complex<double> upperCp = sampleCp[lowerIndex + 1];
            const std::complex<double> lowerCs = sampleCs[lowerIndex];
            const std::complex<double> upperCs = sampleCs[lowerIndex + 1];
            std::complex<double> cp;
            std::complex<double> cs = 0.0;
            if (input.interpolation == AcousticInterpolation::N2Linear)
            {
                const std::complex<double> n2 =
                    (1.0 - fraction) / (lowerCp * lowerCp) +
                    fraction / (upperCp * upperCp);
                cp = 1.0 / std::sqrt(n2);
                if (elastic)
                {
                    const std::complex<double> shearN2 =
                        (1.0 - fraction) / (lowerCs * lowerCs) +
                        fraction / (upperCs * upperCs);
                    cs = 1.0 / std::sqrt(shearN2);
                }
            }
            else
            {
                cp = (1.0 - fraction) * lowerCp + fraction * upperCp;
                if (elastic)
                {
                    cs = (1.0 - fraction) * lowerCs + fraction * upperCs;
                }
            }
            const double rho = (1.0 - fraction) * lower.rho + fraction * upper.rho;

            matrix.depth.push_back(depth);
            matrix.cp.push_back(cp);
            matrix.cs.push_back(cs);
            matrix.rho.push_back(rho);
            if (!elastic)
            {
                matrix.b1.push_back(-2.0 + spacing * spacing * matrix.omega2 / (cp * cp));
                matrix.b2.push_back(0.0);
                matrix.b3.push_back(0.0);
                matrix.b4.push_back(0.0);
                matrix.dynamicRho.push_back(rho);
            }
            else
            {
                const double twoH = 2.0 * spacing;
                const std::complex<double> cp2 = cp * cp;
                const std::complex<double> cs2 = cs * cs;
                matrix.b1.push_back(twoH / (rho * cs2));
                matrix.b2.push_back(twoH / (rho * cp2));
                matrix.b3.push_back(4.0 * twoH * rho * cs2 * (cp2 - cs2) / cp2);
                matrix.b4.push_back(twoH * (cp2 - 2.0 * cs2) / cp2);
                matrix.dynamicRho.push_back(twoH * matrix.omega2 * rho);
            }
        }
        offset += count + 1;
    }
    return matrix;
}
}
