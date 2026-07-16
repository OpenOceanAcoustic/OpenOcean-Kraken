#ifndef OPEN_OCEAN_KRAKENC_PARAMS_H
#define OPEN_OCEAN_KRAKENC_PARAMS_H

#include "OpenOceanKrakencEnvironment.h"

#include <Eigen/Dense>

#include <complex>
#include <string>
#include <vector>

namespace OpenOceanKrakenc
{
constexpr double pi = 3.14159265358979323846;

struct FreqInfo
{
    double freq = 0.0;
    int Nfreq = 0;
    Eigen::VectorXd freqvec;
};

struct MeshParams
{
    int NSets = 5;
    Eigen::VectorXi NV = (Eigen::VectorXi(5) << 1, 2, 4, 8, 16).finished();
};

struct EigenParams
{
    using ComplexScalar = std::complex<double>;

    int M = 0;
    int firstM = 0;
    Eigen::VectorXcd EVMat;
    Eigen::VectorXcd Extrap;
    Eigen::VectorXcd k;
    Eigen::VectorXd VG;
    Eigen::MatrixXcd PsiR;
    Eigen::MatrixXcd PsiS;
    Eigen::MatrixXcd dPsidzR;
    Eigen::MatrixXcd dPsidzS;
    Eigen::VectorXd ModeZ;
    Eigen::MatrixXcd PhiMode;

    void resize(int modalCapacity, int nSets, int nSourceDepths, int nReceiverDepths);
    void setZero();
};

struct TridMtx
{
    using ComplexScalar = std::complex<double>;

    Eigen::VectorXd z;
    Eigen::VectorXcd cp;
    Eigen::VectorXcd cs;
    Eigen::VectorXd rho;
    Eigen::VectorXcd B1;
    Eigen::VectorXcd B2;
    Eigen::VectorXcd B3;
    Eigen::VectorXcd B4;
    Eigen::VectorXi N;
    Eigen::VectorXd h;
    Eigen::VectorXd hV;
    Eigen::VectorXi Loc;
    double cLow = 0.0;
    double cHigh = 0.0;

    void resize(int pointCount, int mediaCount, int nSets);
    void setZero();
};

struct OOKC_parameters
{
    std::string envPath;
    std::string flpPath;
    std::string modPath;
    std::string shdPath;
    std::string Title;
    int totalTasks = 0;
    int NMeshMax = 0;
    int NMediaMax = 0;
    FreqInfo freqinfo;
    Atten_Mode AttenUnit;
    Position Pos;
    Position ModePos;
    bool hasModePos = false;
    int MLimit = 9999;
    int NProf = 1;
    Eigen::VectorXd RProf = Eigen::VectorXd::Zero(1);
    std::vector<ssp::SSPStructure> SSP;
    std::vector<ssp::Range_Independent_Area> sspInput;
    ReflectionCoefInfo ReflectionCoef;
    SrcBmPat SBP;
    bool is_Velocity = false;
    MeshParams mesh;
    double cLow = 0.0;
    double cHigh = 0.0;
    double Rmax = 0.0;
    Source_Mode SourceType = Source_Mode::MODE_R_Point;
    Run_Mode runMode = Run_Mode::MODE_B_Both;
    CoherenceType coherenceType = CoherenceType::Coherent;
    ModeType modeType = ModeType::Adiabatic;
};

struct OOKC_output
{
    std::vector<EigenParams> eigen;
    std::vector<std::complex<float>> pressure;
    std::vector<std::complex<float>> horizontalVelocity;
    std::vector<std::complex<float>> verticalVelocity;
    std::size_t sourceCount = 0;
    std::size_t receiverDepthCount = 0;
    std::size_t rangeCount = 0;

    void clear();
};
}

#endif
