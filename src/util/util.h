#ifndef UTIL_HPP
#define UTIL_HPP
#include "OpenOceanKrakenParams.h"

namespace OpenOceanKraken
{
    namespace Util
    {
        void linspace(double a, double b, int n, Eigen::VectorXd &x);
        void SubTab(Eigen::VectorXd &x, int Nx);
        void Sort(Eigen::VectorXd &x, int N);

        double spacing(double const &x);
        bool isSmallValue(double delta, double ref, bool const_thresh); // 计算是否为小值

        std::complex<float> cpxd2cpxf(const std::complex<double> &a);
        std::complex<double> cpxf2cpxd(const std::complex<float> &a);
    }
}
#endif // UTIL_HPP
