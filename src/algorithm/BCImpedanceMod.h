#ifndef BCImpedanceMod_H
#define BCImpedanceMod_H

#include "kkc_params.h"
#include "RefCoef.h"

// 计算边界条件阻抗
void BCImpedance(const size_t& iprof, const double& x, const bool& isTop, complex<double> &f, complex<double> &g,
                 int &iPower, const bool& isComplex, TridMtx &trid, const parameters& params, 
                 int& modeCount);
void ElasticUP(const double& x, VectorXd &yV, int &iPower, const int& Medium, TridMtx &trid, const parameters& params);
void ElasticDN(const double x, VectorXd &yV, int &iPower, const int& Medium, TridMtx &trid, const parameters& params);

#endif // BCImpedanceMod_H