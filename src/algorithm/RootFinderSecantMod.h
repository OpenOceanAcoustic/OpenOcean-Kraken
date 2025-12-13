#ifndef ROOT_FINDER_SECANT_MOD_H
#define ROOT_FINDER_SECANT_MOD_H

#include "kkc_params.h"

// 定义函数类型，对应Fortran中的FUNCT子例程
typedef void (*RealFunctType)(int& iset, size_t iprof, int &mode, double& x, double &Delta, int &iPower, KrakenMatrix &kramtrx,
           parameters& params, VectorXd &EVMat, const int& firstM, bool coutmodes, int &modeCount);

           // 复数形式待补充
typedef void (*ComplexFunctType)(int& iset, size_t iprof, int &mode, complex<double>& x, complex<double> &Delta, int &iPower, KrakenMatrix &kramtrx,
           parameters& params, VectorXd &EVMat, const int& firstM, bool coutmodes, int &modeCount);



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
void ZSecantX(double &x2, const double Tolerance, int &Iteration, const int MaxIteration,
    int& iset, size_t iprof, int &mode, double &Delta, int &iPower, KrakenMatrix &kramtrx,
           parameters& params, VectorXd &EVMat, const int& firstM, bool coutmodes, int &modeCount,
              std::string &ErrorMessage, RealFunctType Funct);

/**
 * 使用割线法查找复数函数的根
 *
 * @param x2 输入输出参数，初始猜测值，返回找到的根
 * @param Tolerance 输入参数，根的误差界限
 * @param Iteration 输出参数，使用的迭代次数
 * @param MaxIteration 输入参数，最大允许的迭代次数
 * @param ErrorMessage 输出参数，错误信息
 * @param Funct 输入参数，计算函数值的函数指针
 */
void ZSecantCX(std::complex<double> &x2, const double Tolerance, int &Iteration, const int MaxIteration,
    int& iset, size_t iprof, int &mode, double &Delta, int &iPower, KrakenMatrix &kramtrx,
           parameters& params, VectorXd &EVMat, const int& firstM, bool coutmodes, int &modeCount,
               std::string &ErrorMessage, ComplexFunctType Funct);

#endif // ROOT_FINDER_SECANT_MOD_H