#ifndef REFCOEF_H
#define REFCOEF_H

#include "OpenOceanKrakenParams.h"

namespace OpenOceanKraken
{
    void InterpolateReflectionCoefficient(ReflectionCoef &RInt, const Eigen::Matrix<ReflectionCoef, 1, Eigen::Dynamic> &R);
}

#endif // REFCOEF_H