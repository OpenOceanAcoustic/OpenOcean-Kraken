#include "OpenOceanKrakencParams.h"

#include <stdexcept>

namespace OpenOceanKrakenc
{
namespace
{
void requirePositive(int value, const char *name)
{
    if (value <= 0)
    {
        throw std::invalid_argument(std::string(name) + " must be positive");
    }
}

void requireNonNegative(int value, const char *name)
{
    if (value < 0)
    {
        throw std::invalid_argument(std::string(name) + " must not be negative");
    }
}
}

void EigenParams::resize(int modalCapacity, int nSets, int nSourceDepths, int nReceiverDepths)
{
    requirePositive(modalCapacity, "modalCapacity");
    requirePositive(nSets, "nSets");
    requireNonNegative(nSourceDepths, "nSourceDepths");
    requireNonNegative(nReceiverDepths, "nReceiverDepths");

    M = 0;
    firstM = modalCapacity;
    EVMat = Eigen::VectorXcd::Zero(modalCapacity * nSets);
    Extrap = Eigen::VectorXcd::Zero(modalCapacity * nSets);
    k = Eigen::VectorXcd::Zero(modalCapacity);
    VG = Eigen::VectorXd::Zero(modalCapacity);
    PsiS = Eigen::MatrixXcd::Zero(modalCapacity, nSourceDepths);
    PsiR = Eigen::MatrixXcd::Zero(modalCapacity, nReceiverDepths);
    dPsidzS = Eigen::MatrixXcd::Zero(modalCapacity, nSourceDepths);
    dPsidzR = Eigen::MatrixXcd::Zero(modalCapacity, nReceiverDepths);
    ModeZ.resize(0);
    PhiMode.resize(0, 0);
}

void EigenParams::setZero()
{
    M = 0;
    EVMat.setZero();
    Extrap.setZero();
    k.setZero();
    VG.setZero();
    PsiS.setZero();
    PsiR.setZero();
    dPsidzS.setZero();
    dPsidzR.setZero();
    ModeZ.setZero();
    PhiMode.setZero();
}

void TridMtx::resize(int pointCount, int mediaCount, int nSets)
{
    requirePositive(pointCount, "pointCount");
    requirePositive(mediaCount, "mediaCount");
    requirePositive(nSets, "nSets");

    z = Eigen::VectorXd::Zero(pointCount);
    cp = Eigen::VectorXcd::Zero(pointCount);
    cs = Eigen::VectorXcd::Zero(pointCount);
    rho = Eigen::VectorXd::Zero(pointCount);
    B1 = Eigen::VectorXcd::Zero(pointCount);
    B2 = Eigen::VectorXcd::Zero(pointCount);
    B3 = Eigen::VectorXcd::Zero(pointCount);
    B4 = Eigen::VectorXcd::Zero(pointCount);
    N = Eigen::VectorXi::Zero(mediaCount);
    h = Eigen::VectorXd::Zero(mediaCount);
    hV = Eigen::VectorXd::Zero(nSets);
    Loc = Eigen::VectorXi::Zero(mediaCount);
    cLow = 0.0;
    cHigh = 0.0;
}

void TridMtx::setZero()
{
    z.setZero();
    cp.setZero();
    cs.setZero();
    rho.setZero();
    B1.setZero();
    B2.setZero();
    B3.setZero();
    B4.setZero();
    N.setZero();
    h.setZero();
    hV.setZero();
    Loc.setZero();
    cLow = 0.0;
    cHigh = 0.0;
}

void OOKC_output::clear()
{
    eigen.clear();
    pressure.clear();
    horizontalVelocity.clear();
    verticalVelocity.clear();
    sourceCount = 0;
    receiverDepthCount = 0;
    rangeCount = 0;
}
}
