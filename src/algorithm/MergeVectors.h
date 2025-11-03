#ifndef MERGEVECTORS_H
#define MERGEVECTORS_H

#include "kkc_params.h"

void MergeVectors(VectorXd& x, VectorXd& y, VectorXd& z, int& NzTab);
void Weight_dble(VectorXd& x, int Nx, 
                 VectorXd& xTab, int NxTab, 
                 VectorXd& w, VectorXi& Ix);

#endif // MERGEVECTORS_H
