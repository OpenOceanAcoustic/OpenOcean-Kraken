#include "OpenOceanKrakencEnvironment.h"

#include <algorithm>

namespace OpenOceanKrakenc::ssp
{
namespace
{
void setHalfSpace(HSInfo &target, double depth, double alphaR, double alphaI,
                  double betaR, double betaI, double rho)
{
    target.Depth = depth;
    target.alphaR = alphaR;
    target.alphaI = alphaI;
    target.betaR = betaR;
    target.betaI = betaI;
    target.cp = {alphaR, alphaI};
    target.cs = {betaR, betaI};
    target.rho = rho;
    target.BC = BC_Mode::MODE_A_Half_space;
}
}

void Range_Independent_Area::set_Bottom_Line(
    double depth, double alphaR, double alphaI,
    double betaR, double betaI, double rho)
{
    setHalfSpace(HSBot, depth, alphaR, alphaI, betaR, betaI, rho);
}

void Range_Independent_Area::set_Top_Line(
    double depth, double alphaR, double alphaI,
    double betaR, double betaI, double rho)
{
    setHalfSpace(HSTop, depth, alphaR, alphaI, betaR, betaI, rho);
}

FlattenedData Range_Independent_Area::flatten() const
{
    FlattenedData result;
    const Eigen::Index mediaCount = static_cast<Eigen::Index>(layers.size());
    result.NPts.resize(mediaCount);
    result.NMesh.resize(mediaCount);
    result.beta.resize(mediaCount);
    result.ft.resize(mediaCount);
    result.sigma.resize(mediaCount);

    Eigen::Index pointCount = 0;
    for (Eigen::Index medium = 0; medium < mediaCount; ++medium)
    {
        const SSPLayer &layer = layers[static_cast<std::size_t>(medium)];
        result.NPts[medium] = layer.npts;
        result.NMesh[medium] = layer.nmesh;
        result.beta[medium] = layer.beta;
        result.ft[medium] = layer.ft;
        result.sigma[medium] = layer.sigma;
        pointCount += layer.z.size();
    }

    result.z.resize(pointCount);
    result.rho.resize(pointCount);
    result.alphaR.resize(pointCount);
    result.alphaI.resize(pointCount);
    result.betaR.resize(pointCount);
    result.betaI.resize(pointCount);
    Eigen::Index offset = 0;
    for (const SSPLayer &layer : layers)
    {
        const Eigen::Index count = layer.z.size();
        result.z.segment(offset, count) = layer.z;
        result.rho.segment(offset, count) = layer.rho;
        result.alphaR.segment(offset, count) = layer.alphaR;
        result.alphaI.segment(offset, count) = layer.alphaI;
        result.betaR.segment(offset, count) = layer.betaR;
        result.betaI.segment(offset, count) = layer.betaI;
        offset += count;
    }
    return result;
}
}
