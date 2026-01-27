#pragma once

#include "OpenOceanKrakenParams.h"


namespace OpenOceanKraken {
void CSpline(Eigen::VectorXd &TAU, Eigen::MatrixXcd &C, int N, int IBCBEG, int IBCEND, int NDIM);
void VSpline(Eigen::VectorXd &TAU, Eigen::VectorXcd &C, int M, int MDIM, Eigen::VectorXcd &F, int N);
std::complex<double> spline(const std::complex<double> *C, double H);


void SplineALL(Eigen::MatrixXcd &C, int &iSegz, double &H, std::complex<double> &F, std::complex<double> &FX, std::complex<double> &FXX);

};