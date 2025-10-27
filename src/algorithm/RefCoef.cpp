#include "RefCoef.h"

void InterpolateReflectionCoefficient(ReflectionCoef &RInt, Matrix<ReflectionCoef, 1, Dynamic> &R)
{

    int iLeft = 1;
    int iRight = R.size();
    double thetaIntr = RInt.theta;

    if (thetaIntr < R(iLeft - 1).theta)
    {
        RInt.R = 0.0;
        RInt.phi = 0.0;
    }
    else if (thetaIntr > R(iRight - 1).theta)
    {
        RInt.R = 0.0;
        RInt.phi = 0.0;
    }
    else
    {
        while (iLeft != iRight - 1)
        {
            int iMid = (iLeft + iRight) / 2;
            if (R(iMid - 1).theta > thetaIntr)
            {
                iRight = iMid;
            }
            else
            {
                iLeft = iMid;
            }
        }

        double alpha = (RInt.theta - R(iLeft - 1).theta) / (R(iRight - 1).theta - R(iLeft - 1).theta);
        RInt.R = (1 - alpha) * R(iLeft - 1).R + alpha * R(iRight - 1).R;
        RInt.phi = (1 - alpha) * R(iLeft - 1).phi + alpha * R(iRight - 1).phi;
    }
}
