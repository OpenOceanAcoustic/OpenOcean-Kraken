#include "RootFinderSecantMod.h"

void ZSecantX(double &x2, const double &Tolerance, int &Iteration, const int &MaxIteration,
              const int &iset, const size_t iprof, int &mode, double &Delta, int &iPower, TridMtx &trid,
              const parameters &params, VectorXd &EVMat, const int &firstM, const bool &isCountMode, int &modeCount,
              std::string &ErrorMessage, RealFunctType Funct)
{
    int iPower0, iPower1;

    // 检查容差是否为正数
    if (Tolerance <= 0.0)
    {
        ErrorMessage = "Non-positive tolerance specified";
        return; // 在C++中使用return代替STOP
    }

    double x0, x1, shift, f0, f1, cNum, cDen;

    // 初始化第一个点
    x1 = x2 + 10.0 * Tolerance;
    Funct(iset, iprof, mode, x1, f1, iPower1, trid, params, EVMat, firstM, isCountMode, modeCount);

    // 主迭代循环
    for (Iteration = 1; Iteration <= MaxIteration; ++Iteration)
    {
        // 更新点
        x0 = x1;
        f0 = f1;
        iPower0 = iPower1;
        x1 = x2;

        // 计算当前点的函数值
        Funct(iset, iprof, mode, x1, f1, iPower1, trid, params, EVMat, firstM, isCountMode, modeCount);

        // 计算割线法的位移，避免溢出
        cNum = f1 * (x1 - x0);
        cDen = f1 - f0 * std::pow(10.0, iPower0 - iPower1);

        if (std::abs(cNum) >= std::abs(cDen * x1))
        {
            shift = 0.1 * Tolerance;
        }
        else
        {
            shift = cNum / cDen;
        }

        // 更新下一个近似值
        x2 = x1 - shift;

        // 收敛检查
        if (std::abs(x2 - x1) + std::abs(x2 - x0) < Tolerance)
        {
            return; // 收敛成功
        }
    }

    // 迭代失败
    ErrorMessage = "Failure to converge in RootFinderSecant";
}

void ZSecantCX(std::complex<double> &x2, const double Tolerance, int &Iteration, const int MaxIteration,
               const int &iset, const size_t &iprof, const int &mode, double &Delta, int &iPower, TridMtx &trid,
               const parameters &params, VectorXd &EVMat, const int &firstM, const bool &isCountMode, int &modeCount,
               std::string &ErrorMessage, ComplexFunctType Funct)
{

    int iPower0, iPower1;
    // 检查容差是否为正数
    if (Tolerance <= 0.0)
    {
        ErrorMessage = "Non-positive tolerance specified";
        return; // 在C++中使用return代替STOP
    }

    std::complex<double> x0, x1, shift, f0, f1, cNum, cDen;

    // 初始化第一个点
    x1 = x2 + 100.0 * Tolerance;
    Funct(iset, iprof, mode, x1, f1, iPower1, trid, params, EVMat, firstM, isCountMode, modeCount);

    // 主迭代循环
    for (Iteration = 1; Iteration <= MaxIteration; ++Iteration)
    {
        // 更新点
        x0 = x1;
        f0 = f1;
        iPower0 = iPower1;
        x1 = x2;

        // 计算当前点的函数值
        Funct(iset, iprof, mode, x1, f1, iPower1, trid, params, EVMat, firstM, isCountMode, modeCount);

        // 计算割线法的位移，避免溢出
        cNum = f1 * (x1 - x0);
        cDen = f1 - f0 * std::pow(10.0, iPower0 - iPower1);

        if (std::abs(cNum) >= std::abs(cDen * x1))
        {
            shift = 0.1 * Tolerance;
        }
        else
        {
            shift = cNum / cDen;
        }

        // 更新下一个近似值
        x2 = x1 - shift;

        // 收敛检查
        if (std::abs(x2 - x1) + std::abs(x2 - x0) < Tolerance)
        {
            return; // 收敛成功
        }
    }

    // 迭代失败
    ErrorMessage = "Failure to converge in RootFinderSecant";
}