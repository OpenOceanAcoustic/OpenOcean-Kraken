#ifndef SSPMOD_H
#define SSPMOD_H
#include "kkc_params.h"
#include "AttenMod.h"
#include "splinec.h"
#include "pchipMod.h"
void EvaluateSSP(SSPStructure &SSP, SSP_Mode &ssptype);
void n2Linear(SSPStructure &SSP);
void cLinear(SSPStructure &SSP);
void cPCHIP(SSPStructure &SSP);
void cCubic(SSPStructure &SSP);
void Analytic(VectorXcd &cp, VectorXcd &cs, VectorXd &rho, int Medium, int N1);
void UpdateSSPLoss(double &freq, double &freq0, int &NMedia, SSP_Mode &SSPType, Atten_Mode &AttenUnit, SSPStructure *SSPList);
void UpdateHSLoss(double& freq, double& freq0, int& Medium, Atten_Mode &AttenUnit, HSInfo &HSTop, HSInfo &HSBot);
#endif