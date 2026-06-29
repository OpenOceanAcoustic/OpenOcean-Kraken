#ifndef ROOT_FINDER_SECANT_MOD_H
#define ROOT_FINDER_SECANT_MOD_H

#include "OpenOceanKrakenParams.h"
#include "BCImpedanceMod.h"

namespace OpenOceanKraken
{
    /**
     * 使用割线法查找实数函数的根
     *
     * @param x2 输入输出参数，初始猜测值，返回找到的根
     * @param Tolerance 输入参数，根的误差界限
     * @param Iteration 输出参数，使用的迭代次数
     * @param MaxIteration 输入参数，最大允许的迭代次数
     * @param ErrorMessage 输出参数，错误信息
     * @param Funct 输入参数，计算函数值的函数指针
     */
    void ZSecantX(double &x2, const double &Tolerance, int &Iteration, const int &MaxIteration,
                  const int &iset, const size_t iprof, int &mode, double &Delta, int &iPower, TridMtx &trid,
                  const OOK_parameters &params, Eigen::VectorXd &EVMat, const int &firstM, const bool &isCountMode, int &modeCount);

    /**
     * 使用Brent方法在给定区间[A, B]内查找函数的零点X
     *
     * @param x 输出参数，找到的零点
     * @param a 输入输出参数，区间左端点
     * @param b 输入输出参数，区间右端点
     * @param t 输入参数，正的容差
     * @param errorMessage 输出参数，错误信息
     * @param funct 输入参数，计算函数值的函数指针
     *              函数形式：void funct(double x, double &g, int &ipow)，其中g * 10^ipow给出函数值
     */

    bool ZBRENTX(double &x, double &a, double &b, const double &t,
                 const int &iset, const size_t &iprof, const int &mode, double &Delta, int &iPower, TridMtx &trid,
                 const OOK_parameters &params, Eigen::VectorXd &EVMat, const int &firstM, const bool &isCountMode, int &modeCount);

    //工具函数
    void AcousticLayers(const size_t &iprof, double x, double &f, double &g, int &iPower, TridMtx &trid, const OOK_parameters &params, const bool &isCountMode, int &modeCount);
    void FUNCT(const int &iset, const size_t &iprof, const int &mode, double &x, double &Delta, int &iPower, TridMtx &trid,
               const OOK_parameters &params, Eigen::VectorXd &EVMat, const int &firstM, const bool &isCountMode, int &modeCount);
    void Bisection(const int &iset, const size_t &iprof, int &mode, double xMin, double xMax, Eigen::VectorXd &xL, Eigen::VectorXd &xR,
                   TridMtx &trid, const OOK_parameters &params, EigenParams &eigen);

}
#endif // ROOT_FINDER_SECANT_MOD_H
