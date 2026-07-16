#ifndef OPEN_OCEAN_KRAKENC_ENVIRONMENT_H
#define OPEN_OCEAN_KRAKENC_ENVIRONMENT_H

#include "OpenOceanKrakencEnums.h"

#include <Eigen/Dense>

#include <cassert>
#include <complex>
#include <cstddef>
#include <string>
#include <vector>

namespace OpenOceanKrakenc
{
struct HSInfo
{
    double alphaR = 0.0;
    double alphaI = 0.0;
    double betaR = 0.0;
    double betaI = 0.0;
    double beta = 0.0;
    double ft = 0.0;
    std::complex<double> cp{0.0, 0.0};
    std::complex<double> cs{0.0, 0.0};
    double rho = 0.0;
    double Depth = 0.0;
    BC_Mode BC = BC_Mode::MODE_V_Vacuum;
};

namespace ssp
{
struct SSPLayer
{
    int npts = 0;
    int nmesh = 0;
    double beta = 0.0;
    double ft = 0.0;
    double sigma = 0.0;
    Media_Mode Material = Media_Mode::MODE_A_Acoustic;
    Eigen::VectorXd z;
    Eigen::VectorXd rho;
    Eigen::VectorXd alphaR;
    Eigen::VectorXd alphaI;
    Eigen::VectorXd betaR;
    Eigen::VectorXd betaI;

    bool empty() const noexcept { return npts == 0 && z.size() == 0; }
    bool is_valid() const noexcept
    {
        const Eigen::Index count = z.size();
        return count >= 2 && npts == count && rho.size() == count &&
               alphaR.size() == count && alphaI.size() == count &&
               betaR.size() == count && betaI.size() == count;
    }
};

struct FlattenedData
{
    Eigen::VectorXi NPts;
    Eigen::VectorXi NMesh;
    Eigen::VectorXd beta;
    Eigen::VectorXd ft;
    Eigen::VectorXd sigma;
    Eigen::VectorXd z;
    Eigen::VectorXd rho;
    Eigen::VectorXd alphaR;
    Eigen::VectorXd alphaI;
    Eigen::VectorXd betaR;
    Eigen::VectorXd betaI;
};

struct SSPStructure
{
    int NMedia = 0;
    int FirstAcoustic = -1;
    int LastAcoustic = -1;
    Eigen::VectorXd beta;
    Eigen::VectorXd ft;
    Eigen::VectorXd sigma;
    std::vector<Media_Mode> Material;
    SSP_Mode SSPType = SSP_Mode::MODE_C_cLinear;
    Eigen::VectorXd depth;
    Eigen::VectorXi NPts;
    Eigen::VectorXi offset;
    Eigen::VectorXi interp_offset;
    Eigen::VectorXi NMesh;
    Eigen::VectorXd z;
    Eigen::VectorXd alphaR;
    Eigen::VectorXd alphaI;
    Eigen::VectorXd betaR;
    Eigen::VectorXd betaI;
    Eigen::VectorXd rho;
    Eigen::VectorXcd cp;
    Eigen::VectorXcd cs;
    int row_size = 4;
    Eigen::MatrixXcd cspline;
    Eigen::MatrixXcd cCoef;
    Eigen::MatrixXcd csWork;
    Eigen::MatrixXcd cpCoef;
    Eigen::MatrixXcd csCoef;
    Eigen::MatrixXcd rhoCoef;
    Eigen::MatrixXcd cpSpline;
    Eigen::MatrixXcd csSpline;
    Eigen::MatrixXcd rhoSpline;
    HSInfo HSTop;
    HSInfo HSBot;

    int get_media_start(int medium) const
    {
        assert(medium >= 0 && medium < NMedia);
        return offset[medium];
    }
    int get_media_end(int medium) const
    {
        return get_media_start(medium) + NPts[medium];
    }
    int get_media_size(int medium) const
    {
        assert(medium >= 0 && medium < NMedia);
        return NPts[medium];
    }
    int get_media_Nmesh(int medium) const
    {
        assert(medium >= 0 && medium < NMedia);
        return NMesh[medium];
    }
    int get_global_interp_offset(int medium) const
    {
        assert(medium >= 0 && medium < NMedia);
        return interp_offset[medium];
    }
};

struct Range_Independent_Area
{
    std::string Title;
    SSP_Mode SSPType = SSP_Mode::MODE_C_cLinear;
    std::vector<SSPLayer> layers;
    HSInfo HSTop;
    HSInfo HSBot;
    double Range = 0.0;
    bool enableRootRestarts = false;

    void addLayer(const SSPLayer &layer) { layers.push_back(layer); }
    void insertLayer(std::size_t index, const SSPLayer &layer)
    {
        assert(index <= layers.size());
        layers.insert(layers.begin() + static_cast<std::ptrdiff_t>(index), layer);
    }
    void removeLayer(std::size_t index)
    {
        assert(index < layers.size());
        layers.erase(layers.begin() + static_cast<std::ptrdiff_t>(index));
    }
    void clearLayer() { layers.clear(); }
    int NMedia() const noexcept { return static_cast<int>(layers.size()); }

    void set_Bottom_Line(double depth, double alphaR, double alphaI,
                         double betaR, double betaI, double rho);
    void set_Top_Line(double depth, double alphaR, double alphaI,
                      double betaR, double betaI, double rho);
    void set_Bottom_type(BC_Mode type) { HSBot.BC = type; }
    void set_Top_type(BC_Mode type) { HSTop.BC = type; }
    FlattenedData flatten() const;
};
}

struct Position
{
    int NSz = 0;
    int NRz = 0;
    int NRr = 0;
    int NRo = 0;
    int NRz_per_range = 0;
    double Delta_r = 0.0;
    Eigen::VectorXd Sz;
    Eigen::VectorXd Rr;
    Eigen::VectorXd Rz;
    Eigen::VectorXd Ro;
    bool is_Linspace_Rr = false;
    bool is_Linspace_Rz = false;
    bool is_Linspace_Sz = false;
    bool is_Linspace_Ro = false;
    Grid_Mode GridType = Grid_Mode::MODE_R_Rectangular;
};

struct ReflectionCoef
{
    double theta = 0.0;
    double R = 0.0;
    double phi = 0.0;
};

struct InternalReflectionCoefInfo
{
    double freq = 0.0;
    Eigen::VectorXd xTab;
    Eigen::VectorXcd fTab;
    Eigen::VectorXcd gTab;
    Eigen::VectorXi iTab;
    bool isSet = false;
};

struct ReflectionCoefInfo
{
    Eigen::Matrix<ReflectionCoef, 1, Eigen::Dynamic> RTop;
    Eigen::Matrix<ReflectionCoef, 1, Eigen::Dynamic> RBot;
    InternalReflectionCoefInfo IRC;
    bool isDeg = false;
};

struct SrcBmPat
{
    int NSBPPts = 0;
    Eigen::VectorXd theta;
    Eigen::VectorXd pat;
    bool isSet = false;
};

struct BdryType
{
    HSInfo Top;
    HSInfo Bot;
};
}

#endif
