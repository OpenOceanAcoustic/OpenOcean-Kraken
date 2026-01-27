#ifndef SSPMOD_H
#define SSPMOD_H
#include "OpenOceanKrakenParams.h"
#include "AttenMod.h"
#include "splinec.h"
#include "pchipMod.h"
namespace OpenOceanKraken
{
    namespace ssp
    {
        void EvaluateSSP(TridMtx &trid, SSPStructure &SSP, SSP_Mode &ssptype, int iMedium);
        void n2Linear(TridMtx &trid, SSPStructure &SSP, int iMedium);
        void cLinear(TridMtx &trid, SSPStructure &SSP, int iMedium);
        void cPCHIP(TridMtx &trid, SSPStructure &SSP, int iMedium);
        void cCubic(TridMtx &trid, SSPStructure &SSP, int iMedium);
        void Analytic(Eigen::VectorXcd &cp, Eigen::VectorXcd &cs, Eigen::VectorXd &rho, int iMedium, int N1);
        void UpdateSSPLoss(double freq, double freq0,
                           int NMedia, SSP_Mode SSPType, Atten_Mode AttenUnit,
                           SSPStructure &ssp);
        void UpdateHSLoss(double &freq, double &freq0, Atten_Mode &AttenUnit, HSInfo &HSTop, HSInfo &HSBot);
    }
}
#endif