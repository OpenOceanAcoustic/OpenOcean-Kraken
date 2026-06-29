#include "RootFinderMod.h"

namespace OpenOceanKraken
{
    void ZSecantX(double &x2, const double &Tolerance, int &Iteration, const int &MaxIteration,
                  const int &iset, const size_t iprof, int &mode, double &Delta, int &iPower, TridMtx &trid,
                  const OOK_parameters &params, Eigen::VectorXd &EVMat, const int &firstM, const bool &isCountMode, int &modeCount)
    {
        int iPower0, iPower1;

        // 检查容差是否为正数
        if (Tolerance <= 0.0)
        {
            // ErrorMessage = "Non-positive tolerance specified";
            return; // 在C++中使用return代替STOP
        }

        double x0, x1, shift, f0, f1, cNum, cDen;

        // 初始化第一个点
        x1 = x2 + 10.0 * Tolerance;
        FUNCT(iset, iprof, mode, x1, f1, iPower1, trid, params, EVMat, firstM, isCountMode, modeCount);

        // 主迭代循环
        for (Iteration = 1; Iteration <= MaxIteration; ++Iteration)
        {
            // 更新点
            x0 = x1;
            f0 = f1;
            iPower0 = iPower1;
            x1 = x2;

            // 计算当前点的函数值
            FUNCT(iset, iprof, mode, x1, f1, iPower1, trid, params, EVMat, firstM, isCountMode, modeCount);

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

        // // 迭代失败
        // ErrorMessage = "Failure to converge in RootFinderSecant";
    }

    bool ZBRENTX(double &x, double &a, double &b, const double &t,
                 const int &iset, const size_t &iprof, const int &mode, double &Delta, int &iPower, TridMtx &trid,
                 const OOK_parameters &params, Eigen::VectorXd &EVMat, const int &firstM, const bool &isCountMode, int &modeCount)
    {
        int iExpA, iExpB, iExpC;
        const double MACHEP = 1.0E-16;
        const double TEN = 10.0;
        const int MaxIteration = 100;

        double fa, fb, fc, F1, F2, C, D, E, M, P, Q, R, S, TOL;

        // 计算区间端点的函数值
        FUNCT(iset, iprof, mode, a, fa, iExpA, trid, params, EVMat, firstM, isCountMode, modeCount);
        FUNCT(iset, iprof, mode, b, fb, iExpB, trid, params, EVMat, firstM, isCountMode, modeCount);

        // 检查区间端点函数值是否异号
        if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0))
        {
            // errorMessage = "Function sign is the same at the interval endpoints";
            x = a;
            return false;
        }

        // 内部根初始化
        C = a;
        fc = fa;
        iExpC = iExpA;
        E = b - a;
        D = E;
        // 调整值，确保b是当前最佳近似点
        if (iExpA < iExpB)
        {
            F1 = fc * std::pow(TEN, iExpC - iExpB);
            F2 = fb;
        }
        else
        {
            F1 = fc;
            F2 = fb * std::pow(TEN, iExpB - iExpC);
        }

        for (int iteration = 0; iteration < MaxIteration; ++iteration)
        {

            if (std::abs(F1) < std::abs(F2))
            {
                // 交换a, b, c和相应的函数值
                a = b;
                b = C;
                C = a;

                fa = fb;
                iExpA = iExpB;

                fb = fc;
                iExpB = iExpC;

                fc = fa;
                iExpC = iExpA;
            }

            // 计算容差和中点
            TOL = 2.0 * MACHEP * std::abs(b) + t;
            M = 0.5 * (C - b);

            // 检查是否收敛
            if (std::abs(M) <= TOL || fb == 0.0)
            {
                x = b;
                return true;
            }

            // 检查是否需要强制二分法
            if (iExpA < iExpB)
            {
                F1 = fa * std::pow(TEN, iExpA - iExpB);
                F2 = fb;
            }
            else
            {
                F1 = fa;
                F2 = fb * std::pow(TEN, iExpB - iExpA);
            }

            if (std::abs(E) < TOL || std::abs(F1) <= std::abs(F2))
            {
                // 使用二分法
                E = M;
                D = E;
            }
            else
            {
                // 计算插值步长
                S = fb / fa * std::pow(TEN, iExpB - iExpA);

                if (a == C)
                {
                    // 线性插值
                    P = 2.0 * M * S;
                    Q = 1.0 - S;
                }
                else
                {
                    // 逆二次插值
                    Q = fa / fc * std::pow(TEN, iExpA - iExpC);
                    R = fb / fc * std::pow(TEN, iExpB - iExpC);
                    P = S * (2.0 * M * Q * (Q - R) - (b - a) * (R - 1.0));
                    Q = (Q - 1.0) * (R - 1.0) * (S - 1.0);
                }

                // 保证Q为正
                if (P > 0.0)
                {
                    Q = -Q;
                }
                else
                {
                    P = -P;
                }

                S = E;
                E = D;

                // 检查步长是否可接受
                if (2.0 * P < 3.0 * M * Q - std::abs(TOL * Q) && P < std::abs(0.5 * S * Q))
                {
                    D = P / Q;
                }
                else
                {
                    E = M;
                    D = E;
                }
            }

            // 更新a和b
            a = b;
            fa = fb;
            iExpA = iExpB;

            if (std::abs(D) > TOL)
            {
                b += D;
            }
            else
            {
                if (M > 0.0)
                {
                    b += TOL;
                }
                else
                {
                    b -= TOL;
                }
            }

            // 计算新的b点函数值
            FUNCT(iset, iprof, mode, b, fb, iExpB, trid, params, EVMat, firstM, false, modeCount);

            // 更新c点
            if ((fb > 0.0) == (fc > 0.0))
            {
                C = a;
                fc = fa;
                iExpC = iExpA;
                E = b - a;
                D = E;

                // 调整值，确保b是当前最佳近似点
                if (iExpA < iExpB)
                {
                    F1 = fc * std::pow(TEN, iExpC - iExpB);
                    F2 = fb;
                }
                else
                {
                    F1 = fc;
                    F2 = fb * std::pow(TEN, iExpB - iExpC);
                }
            }
        }

        // 返回找到的零点
        x = b;
        return false;
    }

    void AcousticLayers(const size_t &iprof, double x, double &f, double &g, int &iPower, TridMtx &trid, const OOK_parameters &params, const bool &isCountMode, int &modeCount)
    {
        double p0 = 0.0, p1, p2, h2k2;

        if (params.SSP[iprof].FirstAcoustic == -1)
        {
            return;
        }

        // 遍历声学层
        for (int im = params.SSP[iprof].LastAcoustic; im >= params.SSP[iprof].FirstAcoustic; --im)
        {
            h2k2 = SQ(trid.h(im)) * x;
            int ii = trid.Loc(im) + trid.N(im);
            double rhoMedium = trid.rho(trid.Loc(im)); // 使用每层顶部的密度值

            p1 = -2.0 * g;
            p2 = (trid.B1(ii) - h2k2) * g - 2.0 * trid.h(im) * f * rhoMedium;

            // 向上穿过单层介质
            for (ii = trid.Loc(im) + trid.N(im) - 1; ii >= trid.Loc(im); --ii)
            {
                p0 = p1;
                p1 = p2;
                p2 = (h2k2 - trid.B1(ii)) * p1 - p0;

                if (isCountMode)
                {
                    if (p0 * p1 <= 0.0)
                    {
                        modeCount++;
                    }
                }

                // 必要时进行缩放
                if (std::abs(p2) > BCIRoof)
                {
                    p0 = BCIFloor * p0;
                    p1 = BCIFloor * p1;
                    p2 = BCIFloor * p2;
                    iPower = iPower - BCIiPowerF;
                }
            }

            // 计算f和g
            rhoMedium = trid.rho(trid.Loc(im));
            f = -(p2 - p0) / (2.0 * trid.h(im)) / rhoMedium;
            g = -p1;
        }
    }

    void FUNCT(const int &iset, const size_t &iprof, const int &mode, double &x, double &Delta, int &iPower, TridMtx &trid,
               const OOK_parameters &params, Eigen::VectorXd &EVMat, const int &firstM, const bool &isCountMode, int &modeCount)
    {
        int iPowerBot;
        double f, g;
        std::complex<double> fTop, gTop, fBot, gBot;
        bool isTop = false, isComplex = false;

        modeCount = 0;

        // 调用BCImpedance计算底部阻抗
        BCImpedance(iprof, x, isTop, fTop, gTop, iPower, isComplex, trid,
                    params, modeCount);

        f = std::real(fTop);
        g = std::real(gTop);

        // 穿过声学层
        AcousticLayers(iprof, x, f, g, iPower, trid, params, isCountMode, modeCount);

        isTop = true;
        // 调用BCImpedance计算顶部阻抗
        BCImpedance(iprof, x, isTop, fBot, gBot, iPowerBot, isComplex, trid,
                    params, modeCount);

        Delta = std::real(f * std::real(gBot) - g * std::real(fBot));
        iPower = iPower + iPowerBot;

        if (g * Delta > 0.0)
        {
            modeCount++;
        }

        // 减去之前的根
        if (mode > 0 && params.SSP[iprof].NMedia > params.SSP[iprof].LastAcoustic - params.SSP[iprof].FirstAcoustic + 1)
        {
            for (int j = 0; j < mode; ++j)
            {
                Delta = Delta / (x - EVMat(iset * firstM + j));

                // 必要时进行缩放
                while (std::abs(Delta) < BCIFloor && std::abs(Delta) > 0.0)
                {
                    Delta = BCIRoof * Delta;
                    iPower = iPower - BCIiPowerR;
                }

                while (std::abs(Delta) > BCIRoof)
                {
                    Delta = BCIFloor * Delta;
                    iPower = iPower - BCIiPowerF;
                }
            }
        }
    }

    void Bisection(const int &iset, const size_t &iprof, int &mode, double xMin, double xMax, Eigen::VectorXd &xL, Eigen::VectorXd &xR,
                   TridMtx &trid, const OOK_parameters &params, EigenParams &eigen)
    {
        int j, NZeros, NZer1, iPower, modeCount;
        double x, x1, x2, Delta;

        // 初始化左右边界
        for (int i = 0; i < xL.size(); ++i)
        {
            xL(i) = xMin;
            xR(i) = xMax;
        }
        FUNCT(iset, iprof, mode, xMax, Delta, iPower, trid, params, eigen.EVMat, eigen.firstM, true, NZer1);

        if (eigen.M == 1)
        {
            return; // 如果只寻找一个模态，快速退出
        }

        // 遍历每个本征值
        for (int modeIdx = 0; modeIdx < eigen.M - 1; ++modeIdx)
        {
            if (xL(modeIdx) == xMin)
            {
                x2 = xR(modeIdx);

                // 计算x1的初始值
                x1 = xMin;
                for (int i = modeIdx; i < eigen.M; ++i)
                {
                    if (xL(i) > x1)
                    {
                        x1 = xL(i);
                    }
                }
                x1 = std::max(x1, xMin);

                // 开始二分法
                for (j = 0; j < MaxBisections; ++j)
                {
                    x = x1 + (x2 - x1) / 2;
                    FUNCT(iset, iprof, modeIdx, x, Delta, iPower, trid, params, eigen.EVMat, eigen.firstM, true, modeCount);
                    NZeros = modeCount - NZer1;

                    if (NZeros < modeIdx + 1)
                    {
                        // 零点不够多，这是一个新的右边界
                        x2 = x;
                        xR(modeIdx) = x;
                    }
                    else
                    {
                        // 这是一个新的左边界
                        x1 = x;
                        if (xR(NZeros) >= x)
                        {
                            xR(NZeros) = x;
                        }
                        if (xL(NZeros - 1) <= x)
                        {
                            xL(NZeros - 1) = x;
                        }
                    }

                    // 当我们替换了默认的初始xL值，就完成了
                    if (xL(modeIdx) != xMin)
                    {
                        break;
                    }
                }
            }
        }
    }
}
