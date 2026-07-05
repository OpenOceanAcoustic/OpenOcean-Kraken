#ifndef REFCOEF_H
#define REFCOEF_H

#include "OpenOceanKrakenParams.h"

namespace OpenOceanKraken
{
    void InterpolateReflectionCoefficient(ReflectionCoef &RInt, const Eigen::Matrix<ReflectionCoef, 1, Eigen::Dynamic> &R);
    void InterpolateIRC(const std::complex<double> &x, std::complex<double> &f, std::complex<double> &g,
                        int &iPower, const InternalReflectionCoefInfo &irc);
}

#endif // REFCOEF_H
