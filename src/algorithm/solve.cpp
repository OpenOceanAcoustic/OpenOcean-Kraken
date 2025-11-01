#include "solve.h"

// ERROUT函数声明（假设在其他文件中定义）
void ERROUT(const std::string &routine, const std::string &message)
{

};

// Solve1函数：使用Sturm序列分离本征值以及用Brent求根法求得本征值
void Solve1(int& iset, const int& NSets, EigenParams &eigen, KrakenMatrix &kramtrx, parameters &params)
{
    int iPower = 0, NTotal, NzTab = 0, mode = 0;
    double x, x1, x2, xMin, xMax, Eps, Delta;
    bool iscountm = true;
    std::string ErrorMessage;
    VectorXd xL, xR;
    double omega2 = SQ(2 * pi * params.freqinfo->freq);

    int modeCount = 0; // 对应Fortran的.TRUE.

    // 确定模态数量
    xMin = 1.00001 * omega2 / SQ(params.Chigh);

    FUNCT(iset, mode, xMin, Delta, iPower, kramtrx, params, eigen.EVMat, iscountm, modeCount);
    int M = modeCount;
    kramtrx.modeCount = M;

    // 分配xL和xR的内存
    xL.resize(M + 1);
    xR.resize(M + 1);

    if (iset == 0)
    {
        // 分配EVMat等的内存
        eigen.EVMat.resize(NSets * M);
        eigen.Extrap.resize(NSets * M);
        eigen.k.resize(M);
        eigen.VG.resize(M);
        // 检查内存分配（在C++中vector会自动处理内存分配错误）
    }

    xMax = omega2 / SQ(params.Clow); // 最大波数的平方
    FUNCT(iset, mode, xMax, Delta, iPower, kramtrx, params, eigen.EVMat, iscountm, modeCount);

    M = M - modeCount;
    eigen.M = M;

    if (M == 0)
    {
        // 处理无模态情况
        // LRecordLength = 32; // MOD文件片段长度

        // 打开MODFile并写入头部
        // 在C++中需要使用适当的文件I/O

        // 抛出错误
        ERROUT("KRAKEN", "No modes for given phase speed interval");
    }

    // 计算NTotal
    NTotal = 0;
    for (int i = kramtrx.FirstAcoustic; i <= kramtrx.LastAcoustic; ++i)
    {
        NTotal += kramtrx.N(i);
    }

    if (M > NTotal / 5)
    {
        // 输出警告信息
    }

    Bisection(iset, mode, xMin, xMax, xL, xR, kramtrx, params, eigen); // 初始化上下边界

    // 使用ZBRENT精化每个本征值
    iscountm = false;

    for (int modeIdx = 0; modeIdx < M; ++modeIdx)
    {
        x1 = xL(modeIdx);
        x2 = xR(modeIdx);
        Eps = std::abs(x2) * std::pow(10.0, 2.0 - std::numeric_limits<double>::digits10);
        ZBRENTX(x, x1, x2, Eps, iset, mode, Delta, iPower, kramtrx, params, eigen.EVMat, false, modeCount, 
            ErrorMessage, FUNCT); // Brent求根法

        if (!ErrorMessage.empty())
        {
            // 输出警告信息
        }

        eigen.EVMat(iset * M + modeIdx) = x;
    }
}

// FUNCT函数：计算色散关系
void FUNCT(int& iset, int &mode, double& x, double &Delta, int &iPower, KrakenMatrix &kramtrx,
           parameters& params, VectorXd &EVMat, bool coutmodes, int &modeCount)
{
    int iPowerBot;
    double f, g;
    std::complex<double> fTop, gTop, fBot, gBot;

    modeCount = 0;

    // 调用BCImpedance计算底部阻抗
    BCImpedance(x, false, params.HSBot,
                fTop, gTop, iPower, false, params.NMedia,
                params.freqinfo->freq, kramtrx,
                params.ReflectionCoef.RTop, params.ReflectionCoef.RBot, modeCount);

    f = std::real(fTop);
    g = std::real(gTop);

    // 穿过声学层
    AcousticLayers(x, f, g, iPower, kramtrx, coutmodes, modeCount);

    // 调用BCImpedance计算顶部阻抗
    BCImpedance(x, true, params.HSTop,
                fBot, gBot, iPowerBot, false, params.NMedia,
                params.freqinfo->freq, kramtrx,
                params.ReflectionCoef.RTop, params.ReflectionCoef.RBot, modeCount);

    Delta = std::real(f * std::real(gBot) - g * std::real(fBot));
    iPower = iPower + iPowerBot;

    if (g * Delta > 0.0)
    {
        modeCount++;
    }

    // 减去之前的根
    if (mode > 0 && params.NMedia > kramtrx.LastAcoustic - kramtrx.FirstAcoustic + 1)
    {
        for (int j = 0; j < mode; ++j)
        {
            Delta = Delta / (x - EVMat(iset * mode + j));

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

// AcousticLayers函数：穿过声学层
void AcousticLayers(double x, double &f, double &g, int &iPower, KrakenMatrix &kramtrx, bool &coutmodes, int &modeCount)
{
    double p0 = 0.0, p1, p2, h2k2;

    if (kramtrx.FirstAcoustic == -1)
    {
        return;
    }

    // 遍历声学层
    for (int im = kramtrx.LastAcoustic; im >= kramtrx.FirstAcoustic; --im)
    {
        h2k2 = SQ(kramtrx.h(im)) * x;
        int ii = kramtrx.Loc(im) + kramtrx.N(im);
        double rhoMedium = kramtrx.rho(kramtrx.Loc(im)); // 使用每层顶部的密度值

        p1 = -2.0 * g;
        p2 = (kramtrx.B1(ii) - h2k2) * g - 2.0 * kramtrx.h(im) * f * rhoMedium;

        // 向上穿过单层介质
        for (ii = kramtrx.Loc(im) + kramtrx.N(im) - 1; ii >= kramtrx.Loc(im); --ii)
        {
            p0 = p1;
            p1 = p2;
            p2 = (h2k2 - kramtrx.B1(ii)) * p1 - p0;

            if (coutmodes)
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
        rhoMedium = kramtrx.rho(kramtrx.Loc(im));
        f = -(p2 - p0) / (2.0 * kramtrx.h(im)) / rhoMedium;
        g = -p1;
    }
}

// Bisection函数：返回每个本征值的隔离区间
void Bisection(int& iset, int& mode, double xMin, double xMax, VectorXd &xL, VectorXd &xR,
               KrakenMatrix &kramtrx, parameters &params, EigenParams &eigen)
{
    int j, NZeros, NZer1, iPower, modeCount;
    double x, x1, x2, Delta;


    // 初始化左右边界
    for (int i = 0; i < xL.size(); ++i)
    {
        xL(i) = xMin;
        xR(i) = xMax;
    }
    FUNCT(iset, mode, xMax, Delta, iPower, kramtrx, params, eigen.EVMat, true, NZer1);

    if (kramtrx.modeCount == 1)
    {
        return; // 如果只寻找一个模态，快速退出
    }

    // 遍历每个本征值
    for (int modeIdx = 0; modeIdx < kramtrx.modeCount-1; ++modeIdx)
    {
        if (xL(modeIdx) == xMin)
        {
            x2 = xR(modeIdx);

            // 计算x1的初始值
            x1 = xMin;
            for (int i = modeIdx; i < kramtrx.modeCount ; ++i)
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
                FUNCT(iset, modeIdx, x, Delta, iPower, kramtrx, params, eigen.EVMat, true, modeCount);
                NZeros = modeCount - NZer1;

                if (NZeros < modeIdx)
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