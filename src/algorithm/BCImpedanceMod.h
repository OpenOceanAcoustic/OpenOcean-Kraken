#ifndef BCImpedanceMod_H
#define BCImpedanceMod_H

#include "kkc_params.h"
#include "RefCoef.h"

// 计算边界条件阻抗
void BCImpedance(size_t iprof, const double& x,  bool& isTop, complex<double> &f, complex<double> &g,
                 int &iPower, const bool& isComplex, KrakenMatrix &kramtrx, parameters& params, 
                 int& modeCount);
void ElasticUP(const double& x, VectorXd &yV, int &iPower, const int& Medium, KrakenMatrix &kramtrx, parameters& params);
void ElasticDN(const double x, VectorXd &yV, int &iPower, const int& Medium, KrakenMatrix &kramtrx, parameters& params);

#endif // BCImpedanceMod_H