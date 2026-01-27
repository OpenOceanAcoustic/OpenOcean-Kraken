#ifndef MERGEVECTORS_H
#define MERGEVECTORS_H

#include "OpenOceanKrakenParams.h"

namespace OpenOceanKraken
{
    void MergeVectors(const Eigen::VectorXd &x, const Eigen::VectorXd &y, Eigen::VectorXd &z, int &NzTab, Eigen::VectorXi &Ix, Eigen::VectorXi &Iy);
    void Weight_dble(const Eigen::VectorXd &x, const int Nx,
                     const Eigen::VectorXd &xTab, int NxTab,
                     Eigen::VectorXd &w, Eigen::VectorXi &Ix);
}
#endif // MERGEVECTORS_H
