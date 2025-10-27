#pragma once

#include "kkc_params.h"

void CSpline(VectorXd &TAU, MatrixXcd &C, int N, int IBCBEG, int IBCEND, int NDIM);
void VSpline(VectorXd &TAU, VectorXcd &C, int M, int MDIM, VectorXcd &F, int N);
std::complex<double> spline(const std::complex<double> *C, double H);
std::complex<double> splinex(const std::complex<double> *C, double H);
std::complex<double> splinexx(const std::complex<double> *C, double H);
void SplineALL(MatrixXcd &C, int &iSegz, double &H, std::complex<double> &F, std::complex<double> &FX, std::complex<double> &FXX);