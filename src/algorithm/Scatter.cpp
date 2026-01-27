#include "Scatter.h"

// 对应Fortran的ScatterRoot函数：计算散射相关的复数根
namespace OpenOceanKraken
{
    std::complex<double> ScatterRoot(const std::complex<double> &z)
    {
        if (z.real() >= 0.0)
        {
            // 实部非负时，直接对复数开方
            return std::sqrt(z);
        }
        else
        {
            // 实部负时，返回 -i * sqrt(-z)
            return -I1D * std::sqrt(-z);
        }
    }

    // 对应Fortran的KupIng函数：计算Kuperman-Ingenito公式的散射损失扰动
    std::complex<double> KupIng(double sigma,
                                const std::complex<double> &eta1Sq,
                                double rho1,
                                const std::complex<double> &eta2Sq,
                                double rho2,
                                const std::complex<double> &P,
                                const std::complex<double> &U)
    {
        std::complex<double> result = 0.0;

        //  sigma为0时直接返回0，无散射损失
        if (sigma == 0.0)
        {
            return result;
        }

        // 计算复数根eta1、eta2，以及Del
        std::complex<double> eta1 = ScatterRoot(eta1Sq);
        std::complex<double> eta2 = ScatterRoot(eta2Sq);
        std::complex<double> Del = rho1 * eta2 + rho2 * eta1;

        // Del非零时计算散射扰动
        if (std::abs(Del) > 1e-15)
        { // 避免除以零，用极小值判断
            std::complex<double> A11, A12, A21, A22;

            // 计算A11：对应原Fortran公式
            A11 = 0.5 * (eta1Sq - eta2Sq) - (rho2 * eta1Sq - rho1 * eta2Sq) * (eta1 + eta2) / Del;

            // 计算A12：含复数单位i的项
            A12 = I1D * std::pow(rho2 - rho1, 2) * eta1 * eta2 / Del;

            // 计算A21：含负i的项
            A21 = -I1D * std::pow(rho2 * eta1Sq - rho1 * eta2Sq, 2) / (rho1 * rho2 * Del);

            // 计算A22：对应原Fortran公式
            A22 = 0.5 * (eta1Sq - eta2Sq) + (rho2 - rho1) * eta1 * eta2 * (eta1 + eta2) / Del;

            // 计算最终散射扰动结果（原Fortran的KupIng）
            result = -std::pow(sigma, 2) * (-A21 * P * P + (A11 - A22) * P * U + A12 * U * U);
        }

        return result;
    }
}