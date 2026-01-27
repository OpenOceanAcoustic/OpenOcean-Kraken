#ifndef BCImpedanceMod_H
#define BCImpedanceMod_H

#include "OpenOceanKrakenParams.h"
#include "RefCoef.h"

namespace OpenOceanKraken
{
    // 计算边界条件阻抗
    void BCImpedance(const size_t &iprof, const double &x, const bool &isTop, std::complex<double> &f, std::complex<double> &g,
                     int &iPower, const bool &isComplex, TridMtx &trid, const OOK_parameters &params,
                     int &modeCount);
    void ElasticUP(const double &x, Eigen::VectorXd &yV, int &iPower, const int &Medium, TridMtx &trid, const OOK_parameters &params);
    void ElasticDN(const double x, Eigen::VectorXd &yV, int &iPower, const int &Medium, TridMtx &trid, const OOK_parameters &params);

}

#endif // BCImpedanceMod_H