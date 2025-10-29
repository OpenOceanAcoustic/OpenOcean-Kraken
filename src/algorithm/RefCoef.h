#ifndef REFCOEF_H
#define REFCOEF_H

#include "kkc_params.h"

void InterpolateReflectionCoefficient(ReflectionCoef &RInt, Matrix<ReflectionCoef, 1, Dynamic> &R);

#endif // REFCOEF_H