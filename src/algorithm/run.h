#ifndef RUN_H
#define RUN_H

#include "kkc_params.h"
#include "util.h"
#include "algorithm/sspMod.h"
#include "algorithm/solve.h"
#include "algorithm/field.h"

void EigenVWorker(size_t iprof, parameters &params, TridMtx& trid, kkc_output &output);
void ComputePressure(const size_t& iprof, parameters &params, kkc_output &output);


#endif