#include "solve.h"

// ERROUT函数声明（假设在其他文件中定义）
void ERROUT(const std::string &routine, const std::string &message) {

};

void Solve(int &iset, const int &NSets, EigenParams &eigen, KrakenMatrix &kramtrx, parameters &params, double &Error)
{
    double omega2 = SQ(2 * pi * params.freqinfo->freq);
    int iprof = 0;
    if (iprof > 0 && iset < 2 && params.modeType == ModeType::Couple)
    {
        Solve3(iset, NSets, eigen, kramtrx, params);
    }
    else if ((iset < 2) && (params.NMedia <= params.LastAcoustic - params.FirstAcoustic + 1))
    {
        Solve1(iset, NSets, eigen, kramtrx, params);
    }
    else
    {
        Solve2(iset, NSets, eigen, kramtrx, params);
    }
    eigen.Extrap.resize(NSets * eigen.M);
    int start_idx = (iset - 1) * eigen.M;
    int end_idx = start_idx + eigen.M - 1;
    eigen.Extrap.segment(start_idx, eigen.M) = eigen.EVMat.segment(start_idx, eigen.M);

    // 查找满足条件的最小位置
    // Extrap(1, 1:M)对应索引为0到M-1
    int Min_Loc = -1;
    double threshold = omega2 / SQ(params.Chigh);
    for (int mode = 0; mode < eigen.M; ++mode)
    {
        if (eigen.Extrap(mode) > threshold)
        {                       // Extrap(1, mode+1)对应索引mode
            Min_Loc = mode + 1; // Fortran是1基索引
            break;
        }
    }
    if (Min_Loc != -1)
    {
        eigen.M = Min_Loc;
    }

    // 计算NTotal：N(FirstAcoustic : LastAcoustic)的和
    int NTotal = 0;
    for (int i = params.FirstAcoustic; i <= params.LastAcoustic; ++i)
    {
        NTotal += params.mesh.N(i);
    }
    int NTotal1 = NTotal + 1;

    // 如果是第一个网格，计算特征向量
    if (iset == 0)
    {
        Vector(NTotal, NTotal1);
    }

    // 初始化误差和KEY
    Error = 1.0e10;
    int KEY = 2 * eigen.M / 3;

    if (iset > 0)
    {
        // T1 = Extrap(1, KEY)，对应索引KEY-1
        double T1 = eigen.Extrap(KEY);

        // 理查森外推改善精度
        for (int j = iset - 1; j >= 0; --j)
        {
            for (int mode = 0; mode < eigen.M; ++mode)
            {
                // 计算当前(j, mode)和(j+1, mode)在向量中的索引
                int idx_j = j * eigen.M + mode;
                int idx_j1 = (j + 1) * eigen.M + mode;

                double x1 = SQ(params.mesh.NV[j]); // NV是0基存储
                double x2 = SQ(params.mesh.NV[iset]);
                double F1 = eigen.Extrap(idx_j);
                double F2 = eigen.Extrap(idx_j1);

                // 理查森外推公式
                eigen.Extrap(idx_j) = F2 - (F1 - F2) / (x2 / x1 - 1.0);
            }
        }

        // 计算误差
        double T2 = eigen.Extrap(KEY);
        Error = abs(T2 - T1);
    }
}

// Solve1函数：使用Sturm序列分离本征值以及用Brent求根法求得本征值
void Solve1(int &iset, const int &NSets, EigenParams &eigen, KrakenMatrix &kramtrx, parameters &params)
{
    int iPower = 0, NTotal, NzTab = 0, mode = 0;
    double x, x1, x2, xMin, xMax, Eps, Delta;
    bool iscountm = true;
    std::string ErrorMessage;
    VectorXd xL, xR;
    double omega2 = SQ(2 * pi * params.freqinfo->freq);

    int modeCount = 0;

    // 确定模态数量
    xMin = 1.00001 * omega2 / SQ(params.Chigh);

    FUNCT(iset, mode, xMin, Delta, iPower, kramtrx, params, eigen.EVMat, iscountm, modeCount);
    int M = modeCount;
    eigen.M = M;

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
    for (int i = params.FirstAcoustic; i <= params.LastAcoustic; ++i)
    {
        NTotal += params.mesh.N(i);
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
        ZBRENTX(x, x1, x2, Eps, iset, mode, Delta, iPower, kramtrx, params, eigen.EVMat, iscountm, modeCount,
                ErrorMessage, FUNCT); // Brent求根法

        if (!ErrorMessage.empty())
        {
            // 输出警告信息
        }

        eigen.EVMat(iset * M + modeIdx) = x;
    }
}

void Solve2(int &iset, const int &NSets, EigenParams &eigen, KrakenMatrix &kramtrx, parameters &params)
{
    double omega2 = SQ(2 * pi * params.freqinfo->freq), x1, x2, Tolerance, Delta;
    double x = omega2 / SQ(params.Clow);
    int Iteration, MaxIteration = 2000, iPower = 0, modeCount = 0;
    string ErrorMessage;

    bool iscountm = false;

    VectorXd P;
    P.resize(10);

    for (int mode = 0; mode < eigen.M; mode++)
    {
        x = 1.0001 * x;
        if (iset >= 1)
        {
            for (int i = 0; i < iset; i++)
            {
                P(i) = eigen.EVMat(iset * eigen.M + mode);
            }
            if (iset >= 2)
            {
                for (int ii = 0; ii < iset - 1; ii++)
                {
                    for (int j = 0; j < iset - ii - 1; j++)
                    {
                        x1 = SQ(params.mesh.hV(j));
                        x2 = SQ(params.mesh.hV(j + ii));
                        P(j) = ((SQ(params.mesh.hV(iset)) - x2) * P(j) -
                                (SQ(params.mesh.hV(iset)) - x1) * P(j + 1)) /
                               (x1 - x2);
                    }
                }
                x = P(1);
            }
        }
        Tolerance = abs(x) * kramtrx.B1.size() * pow(10.0, (1.0 - std::numeric_limits<double>::digits10));
        ZSecantX(x, Tolerance, Iteration, MaxIteration, iset, mode, Delta, iPower, kramtrx, params, eigen.EVMat, iscountm, modeCount, ErrorMessage, FUNCT);
    }
}

void Solve3(int &iset, const int &NSets, EigenParams &eigen, KrakenMatrix &kramtrx, parameters &params)
{
}

// FUNCT函数：计算色散关系
void FUNCT(int &iset, int &mode, double &x, double &Delta, int &iPower, KrakenMatrix &kramtrx,
           parameters &params, VectorXd &EVMat, bool coutmodes, int &modeCount)
{
    int iPowerBot;
    double f, g;
    std::complex<double> fTop, gTop, fBot, gBot;
    bool isTop = false, isComplex = false;

    modeCount = 0;

    // 调用BCImpedance计算底部阻抗
    BCImpedance(x, isTop, fTop, gTop, iPower, isComplex, kramtrx,
                params, modeCount);

    f = std::real(fTop);
    g = std::real(gTop);

    // 穿过声学层
    AcousticLayers(x, f, g, iPower, kramtrx, params, coutmodes, modeCount);

    isTop = true;
    // 调用BCImpedance计算顶部阻抗
    BCImpedance(x, isTop, fBot, gBot, iPowerBot, isComplex, kramtrx,
                params, modeCount);

    Delta = std::real(f * std::real(gBot) - g * std::real(fBot));
    iPower = iPower + iPowerBot;

    if (g * Delta > 0.0)
    {
        modeCount++;
    }

    // 减去之前的根
    if (mode > 0 && params.NMedia > params.LastAcoustic - params.FirstAcoustic + 1)
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
void AcousticLayers(double x, double &f, double &g, int &iPower, KrakenMatrix &kramtrx, parameters &params, bool &coutmodes, int &modeCount)
{
    double p0 = 0.0, p1, p2, h2k2;

    if (params.FirstAcoustic == -1)
    {
        return;
    }

    // 遍历声学层
    for (int im = params.LastAcoustic; im >= params.FirstAcoustic; --im)
    {
        h2k2 = SQ(params.mesh.h(im)) * x;
        int ii = params.mesh.Loc(im) + params.mesh.N(im);
        double rhoMedium = kramtrx.rho(params.mesh.Loc(im)); // 使用每层顶部的密度值

        p1 = -2.0 * g;
        p2 = (kramtrx.B1(ii) - h2k2) * g - 2.0 * params.mesh.h(im) * f * rhoMedium;

        // 向上穿过单层介质
        for (ii = params.mesh.Loc(im) + params.mesh.N(im) - 1; ii >= params.mesh.Loc(im); --ii)
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
        rhoMedium = kramtrx.rho(params.mesh.Loc(im));
        f = -(p2 - p0) / (2.0 * params.mesh.h(im)) / rhoMedium;
        g = -p1;
    }
}

// Bisection函数：返回每个本征值的隔离区间
void Bisection(int &iset, int &mode, double xMin, double xMax, VectorXd &xL, VectorXd &xR,
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

void Vector(int &iset, int &mode, double xMin, double xMax, VectorXd &xL, VectorXd &xR,
            KrakenMatrix &kramtrx, parameters &params, EigenParams &eigen,
            int &NTotal, int &NTotal1)
{
    int j = 0;
    double h_rho;
    VectorXd z(NTotal1), Phi(NTotal1), d(NTotal1), e(NTotal1+1), zTab;
    z(0) = params.SSP->z(params.FirstAcoustic);
    for (int im = params.FirstAcoustic; im < params.LastAcoustic; im++)
    {
        h_rho = params.mesh.h(im) * kramtrx.rho(params.mesh.Loc(im));
        for (int ii = 1; j <= params.mesh.N(im); j++)
        {
            e(j + ii) = 1.0 / h_rho;
            z(j + ii) = z(j) + params.mesh.h(im) * ii;
            j = j + params.mesh.N(im);
        }
    }
    e(NTotal1) = 1.0 / h_rho;
    MergeVectors(params.Pos->Sz, params.Pos->Rz, zTab);
    // TODO

}