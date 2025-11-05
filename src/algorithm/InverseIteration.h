#ifndef INVERSEITERATION_H
#define INVERSEITERATION_H

#include "kkc_params.h"

void InverseIterationD(int N, 
                       VectorXd& D, 
                       VectorXd& E, 
                       int& IERR, 
                       VectorXd& PhiVector);
                       
#endif