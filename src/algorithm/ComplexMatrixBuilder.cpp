#include "algorithm/ComplexMatrixBuilder.h"

#include "OpenOceanKrakencParams.h"
#include "algorithm/SspInterpolation.h"

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

    for (std::size_t mediumIndex = 0;
         mediumIndex < input.layers.size(); ++mediumIndex)
    {
        const AcousticLayer &layer = input.layers[mediumIndex];
        PreparedComplexSspLayer prepared(
            layer, input.interpolation,
            attenuationContext(
                input, layer.attenuationPower,
                layer.transitionFrequency),
            mediumIndex + 1);
        const bool elastic = prepared.elastic();

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

        for (int local = 0; local <= count; ++local)
        {
            double depth = layer.topDepth + local * spacing;
            if (local == count)
            {
                depth = layer.bottomDepth;
            }
            const InterpolatedMaterial material =
                local == 0
                    ? prepared.top()
                    : (local == count
                           ? prepared.bottom()
                           : prepared.evaluate(depth));
            const std::complex<double> cp = material.cp;
            const std::complex<double> cs = material.cs;
            const double rho = material.rho;

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
