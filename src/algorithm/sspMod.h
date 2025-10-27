#ifndef SSPMOD_H
#define SSPMOD_H
#include "kkc_params.h"
#include "AttenMod.h"
void EvaluateSSP(SSPStructure &SSP);
void n2Linear(SSPStructure &SSP);
void cLinear(SSPStructure &SSP);
void cPCHIP(SSPStructure &SSP);
void cCubic(SSPStructure &SSP);
void Analytic(VectorXcd &cp, VectorXcd &cs, VectorXd &rho, int Medium, int N1);
void PCHIP(VectorXd &x, VectorXcd &y, MatrixXcd &PolyCoef, MatrixXcd &csWork, int istart, int nlen);
void PCHIP(VectorXd &x, VectorXd &y, MatrixXd &PolyCoef, MatrixXd &csWorkd, int istart, int nlen);
void CSpline(VectorXd &TAU, MatrixXcd &C, int N, int IBCBEG, int IBCEND, int NDIM, int istart);
void CSpline(VectorXd &TAU, MatrixXd &C, int N, int IBCBEG, int IBCEND, int NDIM, int istart);
#endif