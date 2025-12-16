#ifndef SSPMOD_H
#define SSPMOD_H
#include "kkc_params.h"
#include "AttenMod.h"
#include "splinec.h"
#include "pchipMod.h"
void EvaluateSSP(SSPStructure &SSP, SSP_Mode &ssptype,int iMedium);
void n2Linear(SSPStructure &SSP,int iMedium);
void cLinear(SSPStructure &SSP,int iMedium);
void cPCHIP(SSPStructure &SSP,int iMedium);
void cCubic(SSPStructure &SSP,int iMedium);
void Analytic(VectorXcd &cp, VectorXcd &cs, VectorXd &rho, int iMedium, int N1);
void UpdateSSPLoss(double freq, double freq0,
                   int NMedia, SSP_Mode SSPType, Atten_Mode AttenUnit,
                   SSPStructure& ssp);
void UpdateHSLoss(double& freq, double& freq0, Atten_Mode &AttenUnit, HSInfo &HSTop, HSInfo &HSBot);
#endif