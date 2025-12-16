#ifndef FIELD_H
#define FIELD_H

#include "kkc_params.h"

inline size_t GetFieldAddr(int32_t isz, int32_t id, int32_t ir, const Position *Pos);
void Evaluate(EigenParams &eigen, parameters &params,int isz,
     std::complex<float> *uAllSources,
     std::complex<float> *uAllSources_vr,
     std::complex<float> *uAllSources_vz);
void field(EigenParams &eigen, parameters &params, std::complex<float> *uAllSources, int isz);
void export_shd(std::string filename, parameters& params, std::complex<float> *uAllSources);

inline size_t GetFieldAddr(int32_t isz, int32_t id, int32_t ir, const Position *Pos)
{
    return ((size_t)isz * (size_t)Pos->NRz_per_range + (size_t)id) * (size_t)Pos->NRr + (size_t)ir;
}


#endif // FIELD_H