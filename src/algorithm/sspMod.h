#ifndef SSPMOD_H
#define SSPMOD_H
#include "kkc_params.h"
void EvaluateSSP(VectorXcd &cp, VectorXcd &cs, VectorXd &rho_k, SSPStructure &SSP, int Medium, int N1, double freq, string task);
void n2Linear(VectorXcd &cp, VectorXcd &cs, VectorXd &rho_k, SSPStructure &SSP, int Medium, int N1);
void cLinear(VectorXcd &cp, VectorXcd &cs, VectorXd &rho_k, SSPStructure &SSP, int Medium, int N1);
void cPCHIP(VectorXcd &cp, VectorXcd &cs, VectorXd &rho_k, SSPStructure &SSP, int Medium, int N1);
void cCubic(VectorXcd &cp, VectorXcd &cs, VectorXd &rho_k, SSPStructure &SSP, int Medium, int N1);
void Analytic(VectorXcd &cp, VectorXcd &cs, VectorXd &rho, int Medium, int N1);
void PCHIP(VectorXd &x, VectorXcd &y, MatrixXcd &PolyCoef, MatrixXcd &csWork, int istart, int nlen);
void PCHIP(VectorXd &x, VectorXd &y, MatrixXd &PolyCoef, MatrixXd &csWorkd, int istart, int nlen);
std::complex<double> CRCI(double &z, double &c, double &alpha, double &freq, double &freq0,
                          const Atten_Mode &AttenUnit, double &beta, double &fT);
void CSpline(VectorXd &TAU, MatrixXcd &C, int N, int IBCBEG, int IBCEND, int NDIM, int istart);
void CSpline(VectorXd &TAU, MatrixXd &C, int N, int IBCBEG, int IBCEND, int NDIM, int istart);
#endif