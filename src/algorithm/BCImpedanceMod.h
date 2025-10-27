#ifndef BCImpedanceMod_H
#define BCImpedanceMod_H

#include "kkc_params.h"
#include "RefCoef.h"

void ElasticUP(const double x, VectorXd& yV, int& iPower, const int Medium, KrakenMatrix& kramtrx);
void ElasticDN(const double x, VectorXd& yV, int& iPower, const int Medium, KrakenMatrix& kramtrx);

#endif // BCImpedanceMod_H