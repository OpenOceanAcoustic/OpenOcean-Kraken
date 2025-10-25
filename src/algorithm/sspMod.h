#ifndef SSPMOD_H
#define SSPMOD_H
#include "kkc_params.h"
void EvaluateSSP(VectorXcd &cP, VectorXcd &cS, VectorXd &rho_k, SSPStructure &SSP, int Medium, int N1, double freq, string task);
void n2Linear(VectorXcd &cP, VectorXcd &cS, VectorXd &rho_k, SSPStructure &SSP, int Medium, int N1);
void cLinear(VectorXcd &cP, VectorXcd &cS, VectorXd &rho_k, SSPStructure &SSP, int Medium, int N1);
void cPCHIP(VectorXcd &cP, VectorXcd &cS, VectorXd &rho_k, SSPStructure &SSP, int Medium, int N1);
void cCubic(VectorXcd &cP, VectorXcd &cS, VectorXd &rho_k, SSPStructure &SSP, int Medium, int N1);
void Analytic(VectorXcd& cP, VectorXcd& cS, VectorXd& rho, int Medium, int N1);
#endif