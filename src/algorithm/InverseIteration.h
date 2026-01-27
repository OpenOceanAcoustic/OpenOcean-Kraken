#ifndef INVERSEITERATION_H
#define INVERSEITERATION_H

#include "OpenOceanKrakenParams.h"

namespace OpenOceanKraken
{
void InverseIterationD(int N, 
                       Eigen::VectorXd& D, 
                       Eigen::VectorXd& E, 
                       int& IERR, 
                       Eigen::VectorXd& PhiVector);
}
#endif