#ifndef ROOT_FINDER_BRENT_H
#define ROOT_FINDER_BRENT_H

#include "kkc_params.h"

// 定义函数类型，对应Fortran中的FUNCT子例程
typedef void (*FunctType)(const int &iset, const size_t &iprof, const int &mode, double &x, double &Delta, int &iPower, TridMtx &trid,
                          parameters &params, VectorXd &EVMat, const int &firstM, const bool &isCountMode, int &modeCount);

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

void ZBRENTX(double &x, double &a, double &b, const double t,
             const int &iset, const size_t &iprof, const int &mode, double &Delta, int &iPower, TridMtx &trid,
             parameters &params, VectorXd &EVMat, const int &firstM, int &modeCount,
             std::string &errorMessage, FunctType funct);

#endif // ROOT_FINDER_BRENT_H