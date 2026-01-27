#ifndef FIELD_H
#define FIELD_H

#include "OpenOceanKrakenParams.h"

namespace OpenOceanKraken
{
    inline size_t GetFieldAddr(int32_t isz, int32_t id, int32_t ir, const Position *Pos);
    void Evaluate(EigenParams &eigen, const OOK_parameters &params, int isz,int iprof,
                  std::complex<float> *uAllSources,
                  std::complex<float> *uAllSources_vr,
                  std::complex<float> *uAllSources_vz);
    void field(EigenParams &eigen, const OOK_parameters &params, std::complex<float> *uAllSources, int isz);
    
    void getSourceEnv(const ssp::SSPStructure &ssp,const double Sz, double &rho,double &c0);//插值计算声源位置的声速和密度






    inline size_t GetFieldAddr(int32_t isz, int32_t id, int32_t ir, const Position *Pos)
    {
        return (static_cast<size_t>(isz) * static_cast<size_t>(Pos->NRz_per_range) + static_cast<size_t>(id)) * static_cast<size_t>(Pos->NRr) + static_cast<size_t>(ir);
    }

}

#endif // FIELD_H