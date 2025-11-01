#ifndef BCImpedanceMod_H
#define BCImpedanceMod_H

#include "kkc_params.h"
#include "RefCoef.h"

// 计算边界条件阻抗
void BCImpedance(const double x, bool isTop, const HSInfo &HS,
                 std::complex<double> &f, std::complex<double> &g,
                 int &iPower, const bool ComplexFlag, int &NMedia,
                 double &freq, KrakenMatrix &kramtrx,
                 Matrix<ReflectionCoef, 1, Dynamic> &RTop, Matrix<ReflectionCoef, 1, Dynamic> &RBot, int &modeCount);
void ElasticUP(const double x, VectorXd &yV, int &iPower, const int Medium, KrakenMatrix &kramtrx);
void ElasticDN(const double x, VectorXd &yV, int &iPower, const int Medium, KrakenMatrix &kramtrx);

#endif // BCImpedanceMod_H