#include "run.h"

// 计算本征值和本征函数
void EigenVWorker(ThreadPool &threadPool, const int &NumThreads, const size_t &iprof, const parameters &params, TridMtx& trid, kkc_output &output)
{
    double freq = params.freqinfo->freq;
    auto &ssp = params.SSP[iprof];
    auto &HSTop = params.HSTop[iprof];
    auto &HSBot = params.HSBot[iprof];
    auto &eigen = output.eigen[iprof];
    double error;

    for (int iset = 0; iset < params.mesh.NSets; iset++)
    {
        int ntimes = params.mesh.NV[iset];
        TridPreprocess(iset, iprof, params, trid, ntimes);
        SolveEp(threadPool, NumThreads, iset, iprof, params.mesh.NSets, eigen, trid, params, error);
        if (error * params.Rmax < 1.0)
        {
            break;
        }
        else
        {
            if (iset == params.mesh.NSets - 1)
                cout << "Warning in KRAKEN : Too many meshes needed: check convergence" << endl;
        }
    }

    size_t m = 0;
    while (m < eigen.M && eigen.Extrap(m) > SQ(2 * pi * params.freqinfo->freq / trid.cHigh))
    {
        m++;
    }
    eigen.M = m;

    for (int i = 0; i < m; i++)
    {
        eigen.k(i) = sqrt(eigen.Extrap(i) + eigen.k(i));
    }
}

// 计算声场
void FieldWorker(const size_t& iprof, const parameters &params, kkc_output &output)
{
    auto eigen = output.eigen[iprof];
    for (int isz = 0; isz < params.Pos->NSz; isz++)
    {
        Evaluate(eigen, params, isz, output.u_AllSources, output.v_AllSources, output.h_AllSources);
    }
}

// ERROUT函数声明（假设在其他文件中定义）
void ERROUT() {

};

void SolveEp(ThreadPool &threadPool, const int &NumThreads, const int &iset, const size_t &iprof, const int &NSets, EigenParams &eigen, TridMtx &trid, const parameters &params, double &Error)
{
    double omega2 = SQ(2 * pi * params.freqinfo->freq);
    if (iprof > 0 && iset < 2 && params.modeType == ModeType::Couple)
    {
        Solve3(threadPool, NumThreads, iset, iprof, eigen, trid, params);
    }
    else if ((iset < 2) && (params.SSP[iprof].NMedia <= params.SSP[iprof].LastAcoustic - params.SSP[iprof].FirstAcoustic + 1))
    {
        Solve1(threadPool, NumThreads, iset, iprof, NSets, eigen, trid, params);
    }
    else
    {
        Solve2(threadPool, NumThreads, iset, iprof, eigen, trid, params);
    }

    int start_idx = iset * eigen.firstM;
    eigen.Extrap.segment(start_idx, eigen.M) = eigen.EVMat.segment(start_idx, eigen.M);

    // std::cout << "iset: " << iset << " \n"
    //           << eigen.EVMat.segment(start_idx, eigen.M) << std::endl;

    // 查找满足条件的最小位置
    double threshold = omega2 / SQ(trid.cHigh);

    int Min_Loc = 0;

    while (Min_Loc < eigen.M && eigen.Extrap(Min_Loc) > threshold)
    {
        Min_Loc++;
    }
    if (Min_Loc != 0)
    {
        eigen.M = Min_Loc;
    }

    // 计算NTotal：N(FirstAcoustic : LastAcoustic)的和
    int NTotal = 0;
    for (int i = params.SSP[iprof].FirstAcoustic; i <= params.SSP[iprof].LastAcoustic; ++i)
    {
        NTotal += trid.N(i);
    }
    int NTotal1 = NTotal + 1;

    // 如果是第一个网格，计算特征向量
    if (iset == 0)
    {
        VectorSolve(iprof, trid, params, eigen, NTotal, NTotal1);
    }

    // // 打印特征值
    // std::cout << "iset: " << iset << " \n" << std::endl;
    // for (int mode = 0; mode < eigen.M; ++mode)
    // {
    //     // 8位精度输出
    //     std::cout << std::setprecision(8) << eigen.EVMat(iset * eigen.firstM + mode) << std::endl;
    // }

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
                int idx_j = j * eigen.firstM + mode;
                int idx_j1 = (j + 1) * eigen.firstM + mode;

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
void Solve1(ThreadPool &threadPool, const int &NumThreads, const int &iset, const size_t &iprof, const int &NSets, EigenParams &eigen, TridMtx &trid, const parameters &params)
{
    int iPower = 0, NTotal, mode = 0; // NzTab = 0,
    double x, x1, x2, xMin, xMax, Eps, Delta;
    bool isCountMode = true;
    std::string ErrorMessage;
    VectorXd xL, xR;
    double omega2 = SQ(2 * pi * params.freqinfo->freq);

    int modeCount = 0;

    // 确定模态数量
    xMin = 1.00001 * omega2 / SQ(trid.cHigh);

    FUNCT(iset, iprof, mode, xMin, Delta, iPower, trid, params, eigen.EVMat, eigen.firstM, isCountMode, modeCount);
    int M = modeCount;
    eigen.M = M;

    // 分配xL和xR的内存
    xL.resize(M + 1);
    xR.resize(M + 1);

    xMax = omega2 / SQ(trid.cLow); // 最大波数的平方
    FUNCT(iset, iprof, mode, xMax, Delta, iPower, trid, params, eigen.EVMat, eigen.firstM, isCountMode, modeCount);

    M = M - modeCount;
    eigen.M = M;

    if (M == 0)
    {
        // 处理无模态情况
        // LRecordLength = 32; // MOD文件片段长度

        // 打开MODFile并写入头部
        // 在C++中需要使用适当的文件I/O

        // 抛出错误
        ERROUT();
    }

    // 计算NTotal
    NTotal = 0;
    for (int i = params.SSP[iprof].FirstAcoustic; i <= params.SSP[iprof].LastAcoustic; ++i)
    {
        NTotal += trid.N(i);
    }

    if (M > NTotal / 5)
    {
        // 输出警告信息
        std::cout << "Approximate number of modes = " << M << " Warning in KRAKEN - Solve1 : Mesh too coarse to sample the modes adequately" << std::endl;
    }

    Bisection(iset, iprof, mode, xMin, xMax, xL, xR, trid, params, eigen); // 初始化上下边界

    // 使用ZBRENT精化每个本征值
    isCountMode = false;
    std::cout << "NumThreads: " << NumThreads << std::endl;

    // 并行处理每个模态
    for (int threadId = 0; threadId < NumThreads; ++threadId)
    {
        threadPool.enqueue_with_id("RunZBRENTX", [&NumThreads, &params, &eigen, &trid, &xL, &xR, M, iset, iprof, mode, isCountMode, threadId]()
        {
            // 当前线程处理的任务 ID 满足 modeIdx % NumThreads == threadId
            for (int modeIdx = threadId; modeIdx < M; modeIdx += NumThreads)
            {
                double local_x, local_x1, local_x2, local_Eps, local_Delta;
                int local_iPower = 0, local_modeCount = 0;
                std::string local_ErrorMessage;

                local_x1 = xL(modeIdx);
                local_x2 = xR(modeIdx);
                local_Eps = std::abs(local_x2) * std::pow(10.0, 2.0 - std::numeric_limits<double>::digits10);
                
                ZBRENTX(local_x, local_x1, local_x2, local_Eps, iset, iprof, mode, local_Delta, local_iPower, trid, params, eigen.EVMat, eigen.firstM, isCountMode, local_modeCount,
                        local_ErrorMessage); // Brent求根法

                if (!local_ErrorMessage.empty())
                {
                    // 输出警告信息
                }

                eigen.EVMat(iset * eigen.firstM + modeIdx) = local_x;
            }
        });
    }
    
    // 等待所有ZBRENTX任务完成
    threadPool.wait_id("RunZBRENTX");
}

void Solve2(ThreadPool &threadPool, const int &NumThreads, const int &iset, const size_t &iprof, EigenParams &eigen, TridMtx &trid, const parameters &params)
{
    double omega2 = SQ(2 * pi * params.freqinfo->freq), x1, x2, Tolerance, Delta;
    double x = omega2 / SQ(trid.cLow);
    int Iteration, MaxIteration = 2000, iPower = 0, modeCount = 0;
    string ErrorMessage;

    bool isCountMode = false;

    VectorXd P;
    P.resize(10);

    for (int mode = 0; mode < eigen.M; mode++)
    {
        x = 1.00001 * x;
        if (iset >= 1)
        {
            for (int i = 0; i < iset; i++)
            {
                P(i) = eigen.EVMat(i * eigen.firstM + mode);
            }
            if (iset >= 2)
            {
                for (int ii = 0; ii < iset - 1; ii++)
                {
                    for (int j = 0; j < iset - ii - 1; j++)
                    {
                        x1 = SQ(trid.h(j));
                        x2 = SQ(trid.h(j + ii + 1));
                        P(j) = ((SQ(trid.h(iset)) - x2) * P(j) -
                                (SQ(trid.h(iset)) - x1) * P(j + 1)) /
                               (x1 - x2);
                    }
                }
                x = P(0);
            }
        }
        Tolerance = abs(x) * trid.B1.size() * pow(10.0, (1.0 - std::numeric_limits<double>::digits10));
        ZSecantX(x, Tolerance, Iteration, MaxIteration, iset, iprof, mode, Delta, iPower, trid, params, eigen.EVMat, eigen.firstM, isCountMode, modeCount, ErrorMessage);
        eigen.EVMat(iset * eigen.firstM + mode) = x;
        if (omega2 / SQ(trid.cHigh) > x)
        {
            eigen.M = mode;
            return;
        }
    }
}

void Solve3(ThreadPool &threadPool, const int &NumThreads, const int &iset, const size_t &iprof, EigenParams &eigen, TridMtx &trid, const parameters &params)
{
    int IT, MaxIT, iPower = 0, mode = 0; // NzTab = 0,
    double x, xMin, Tolerance, Delta;
    std::string ErrorMessage;

    bool isCountMode = false;

    MaxIT = 500;

    double omega2 = SQ(2 * pi * params.freqinfo->freq);

    int modeCount = 0;

    // 确定模态数量
    xMin = 1.00001 * omega2 / SQ(trid.cHigh);

    FUNCT(iset, iprof, mode, xMin, Delta, iPower, trid, params, eigen.EVMat, eigen.firstM, isCountMode, modeCount);
    int M = modeCount;

    for (int modeIdx = 0; modeIdx < M; ++modeIdx)
    {
        x = eigen.EVMat(iset * eigen.firstM + modeIdx);
        Tolerance = std::abs(x) * std::pow(10.0, 2.0 - std::numeric_limits<double>::digits10);
        ZSecantX(x, Tolerance, IT, MaxIT, iset, iprof, mode, Delta, iPower, trid, params, eigen.EVMat, eigen.firstM, isCountMode, modeCount, ErrorMessage);

        if (!ErrorMessage.empty())
        {
            // 输出警告信息
        }

        eigen.EVMat(iset * eigen.firstM + modeIdx) = x;

        if (omega2 / SQ(trid.cHigh) > x)
        {
            eigen.M = modeIdx; // 调整为当前索引
            return;
        }
    }
}

// 初始化有限差分方程
void TridPreprocess(int &iset, size_t iprof, const parameters &params, TridMtx &trid, int ntimes)
{
    bool ElasticFlag = false;
    int NPoints = 0;

    double Two_h;
    double cp2, cs2;
    double omega2 = SQ(2 * pi * params.freqinfo->freq);

    // 初始化变量
    double cLow = params.cLow;
    double cMin = 1e8;
    double cHigh = params.cHigh;
    trid.Loc[0] = 0; // C++使用0-based索引

    for (int i = 0; i < params.SSP[iprof].NMedia; ++i)
    {
        trid.N(i) = params.SSP[iprof].NMesh(i) * ntimes;
        trid.h(i) = params.SSP[iprof].depth(i) / trid.N(i);
        // 计算当前层的起始位置
        if (i != 0)
        {
            trid.Loc(i) = trid.Loc(i - 1) + trid.N(i - 1) + 1;
        }
        NPoints += trid.N(i);
        if (i == 0)
        {
            trid.hV(iset) = trid.h(i);
        }
    }
    NPoints += params.SSP[iprof].NMedia;

    // 处理每个介质层
    for (int im = 0; im < params.SSP[iprof].NMedia; ++im)
    { 
        int ii = trid.Loc(im); // C++使用0-based索引，不需要+1
        Two_h = 2.0 * trid.h(im);

        EvaluateSSP(trid, params.SSP[iprof], params.SSP[iprof].SSPType, im);

        // 打印rho_int
        // std::cout << "RHO INT" << SSP.rho_int.transpose() << std::endl;

        // 加载有限差分方程的对角线
        if (params.SSP[iprof].Material[im] == Media_Mode::MODE_A_Acoustic)
        { // 声学介质情况

            // 计算当前层的最小声速
            double min_cp = 1e8;
            for (int j = 0; j < trid.N(im) + 1; ++j)
            {
                min_cp = std::min(min_cp, std::real(trid.cp_int(ii + j)));
            }
            cMin = std::min(cMin, min_cp);

            // 计算B1和B1C
            double h2 = SQ(trid.h(im));
            for (int j = 0; j < trid.N(im) + 1; ++j)
            {
                complex<double> val = omega2 / SQ(trid.cp_int(ii + j));
                trid.B1(ii + j) = -2.0 + h2 * std::real(val);
                trid.B1C(ii + j) = std::imag(val);
                trid.rho(ii + j) = trid.rho_int(ii + j);
            }
        }
        else
        { // 弹性介质情况
            ElasticFlag = true;
            Two_h = 2.0 * trid.h(im);

            for (int j = 0; j < trid.N(im) + 1; ++j)
            {
                cMin = std::min(std::real(trid.cs_int[j]), cMin);

                cp2 = SQ(std::real(trid.cp_int(j)));
                cs2 = SQ(std::real(trid.cs_int(j)));

                trid.B1(ii + j) = Two_h / (trid.rho_int(j) * cs2);
                trid.B2(ii + j) = Two_h / (trid.rho_int(j) * cp2);
                trid.B3(ii + j) = 4.0 * Two_h * trid.rho_int(j) * cs2 * (cp2 - cs2) / cp2;
                trid.B4(ii + j) = Two_h * (cp2 - 2.0 * cs2) / cp2;
                trid.rho(ii + j) = Two_h * std::real(omega2) * trid.rho_int(j);
            }
        }
    }

    // // 打印kramtrx内部参数
    // std::cout << "B1 " << kramtrx.B1.size() << kramtrx.B1.transpose() << std::endl;

    // std::cout << "B1C " << kramtrx.B1C.size() << kramtrx.B1C.transpose() << std::endl;

    // std::cout << "rho " << kramtrx.rho.size() << kramtrx.rho.transpose() << std::endl;

    // 处理底部半空间属性
    if (params.HSBot[iprof].BC == BC_Mode::MODE_A_Half_space)
    {
        if (std::real(params.HSBot[iprof].cs) > 0.0)
        { // 弹性底部半空间
            ElasticFlag = true;
            cMin = std::min(cMin, std::real(params.HSBot[iprof].cs));
            cHigh = std::min(cHigh, std::real(params.HSBot[iprof].cs));
        }
        else
        { // 声学底部半空间
            cMin = std::min(cMin, std::real(params.HSBot[iprof].cp));
        }
    }

    // 处理顶部半空间属性
    if (params.HSTop[iprof].BC == BC_Mode::MODE_A_Half_space)
    {
        if (std::real(params.HSTop[iprof].cs) > 0.0)
        { // 弹性顶部半空间
            ElasticFlag = true;
            cMin = std::min(cMin, std::real(params.HSTop[iprof].cs));
            cHigh = std::min(cHigh, std::real(params.HSTop[iprof].cs));
        }
        else
        { // 声学顶部半空间
            cMin = std::min(cMin, std::real(params.HSTop[iprof].cp));
        }
    }

    // 如果存在弹性介质，则减小cMin以考虑Scholte波
    if (ElasticFlag)
    {
        cMin = 0.85 * cMin;
    }
    trid.cLow = std::max(cLow, cMin);
    trid.cHigh = cHigh;
}

// FUNCT函数：计算色散关系
void FUNCT(const int &iset, const size_t &iprof, const int &mode, double &x, double &Delta, int &iPower, TridMtx &trid,
           const parameters &params, VectorXd &EVMat, const int &firstM, const bool &isCountMode, int &modeCount)
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

// AcousticLayers函数：穿过声学层
void AcousticLayers(const size_t &iprof, double x, double &f, double &g, int &iPower, TridMtx &trid, const parameters &params, const bool &isCountMode, int &modeCount)
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

// Bisection函数：返回每个本征值的隔离区间
void Bisection(const int &iset, const size_t &iprof, int &mode, double xMin, double xMax, VectorXd &xL, VectorXd &xR,
               TridMtx &trid, const parameters &params, EigenParams &eigen)
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

void VectorSolve(size_t iprof, TridMtx &trid, const parameters &params, EigenParams &eigen,
                 int &NTotal, int &NTotal1)
{
    int j = 0, NzTab, iPower, modeCount = 0, L, ITP, IErr;
    double h_rho, x, xh2;
    VectorXd z(NTotal1), Psi(NTotal1), dPsidz(NTotal1), d(NTotal1), e(NTotal1 + 1), zTab, WTS, WTR, WTZ;
    VectorXi ISzTab, IRzTab, IZzTab, Ix, Iy;
    VectorXcd PhiTab;
    complex<double> fTop, gTop, fBot, gBot;
    bool isTop = true, isComplex = false;
    z(0) = params.SSP[iprof].z(params.SSP[iprof].FirstAcoustic);
    for (int im = params.SSP[iprof].FirstAcoustic; im <= params.SSP[iprof].LastAcoustic; im++)
    {
        h_rho = trid.h(im) * trid.rho(trid.Loc(im));
        for (int ii = 1; ii <= trid.N(im); ii++)
        {
            e(j + ii) = 1.0 / h_rho;
            z(j + ii) = z(j) + trid.h(im) * ii;
        }
        j = j + trid.N(im);
    }
    e(NTotal1) = 1.0 / h_rho;
    // 为了方便输出.mod，这里还是进行合并
    MergeVectors(params.Pos->Sz, params.Pos->Rz, zTab, NzTab, Ix, Iy);
    MatrixXcd phiZ(eigen.firstM, NzTab);

    WTS.resize(params.Pos->NSz);
    WTR.resize(params.Pos->NRz);
    WTZ.resize(NzTab);

    ISzTab.resize(params.Pos->NSz);
    IRzTab.resize(params.Pos->NRz);
    IZzTab.resize(NzTab);

    Weight_dble(z, NTotal1, params.Pos->Sz, params.Pos->NSz, WTS, ISzTab);
    Weight_dble(z, NTotal1, params.Pos->Rz, params.Pos->NRz, WTR, IRzTab);
    Weight_dble(z, NTotal1, zTab, NzTab, WTZ, IZzTab);

    for (int mode = 0; mode < eigen.M; mode++)
    {
        x = eigen.EVMat(mode);
        BCImpedance(iprof, x, isTop, fTop, gTop, iPower, isComplex, trid,
                    params, modeCount);
        if (gTop == 0.0)
        {
            d(0) = 1.0;
            e(1) = 0.0;
        }
        else
        {
            L = trid.Loc(params.SSP[iprof].FirstAcoustic);
            xh2 = x * SQ(trid.h(params.SSP[iprof].FirstAcoustic));
            h_rho = trid.h(params.SSP[iprof].FirstAcoustic) * trid.rho(L);
            d(0) = (trid.B1(L) - xh2) / h_rho / 2.0 + real(fTop / gTop);
        }

        ITP = NTotal;
        j = 0;
        L = trid.Loc(params.SSP[iprof].FirstAcoustic);
        for (int im = params.SSP[iprof].FirstAcoustic; im <= params.SSP[iprof].LastAcoustic; im++)
        {
            xh2 = x * SQ(trid.h(im));
            h_rho = trid.h(im) * trid.rho(L + 1);
            if (im >= params.SSP[iprof].FirstAcoustic + 1)
            {
                L++;
                d(j) = (d(j) + (trid.B1(L) - xh2) / h_rho) / 2.0;
            }
            for (int ii = 0; ii < trid.N(im); ii++)
            {
                j++;
                L++;
                d(j) = (trid.B1(L) - xh2) / h_rho;
                if (trid.B1(L) - xh2 + 2.0 > 0.0)
                    ITP = min(j, ITP);
            }
        }
        isTop = false;
        BCImpedance(iprof, x, isTop, fBot, gBot, iPower, isComplex, trid,
                    params, modeCount);
        if (gBot == 0.0)
        {
            d(NTotal1 - 1) = 1.0;
            e(NTotal1 - 1) = 0.0;
        }
        else
        {
            d(NTotal1 - 1) = d(NTotal1 - 1) / 2.0 - real(fBot / gBot);
        }
        InverseIterationD(NTotal1, d, e, IErr, Psi);

        if (IErr != 0)
        {
            // 未收敛
            std::cout << "mode = " << mode << std::endl;
            std::cout << "Warning in KRAKEN - InverseIteration: Inverse iteration failed to converge" << std::endl;
            Psi = VectorXd::Zero(NTotal1); // zero out the errant eigenvector
        }
        else
        {
            Normalize(iprof, mode, eigen.firstM, Psi, ITP, NTotal1, x, trid, params, eigen);
        }
        int index = 0;
        for (int im = params.SSP[iprof].FirstAcoustic; im <= params.SSP[iprof].LastAcoustic; ++im)
        {
            double h = trid.h(im);
            int Nn = trid.N(im);
            if (im == params.SSP[iprof].LastAcoustic)
            {
                Nn++;
            }
            for (int ii = 0; ii < Nn; ii++)
            {
                if (index == 0)
                {
                    dPsidz(index) = (Psi(index + 1) - Psi(index)) / h;
                }
                else if (index == NTotal)
                {
                    dPsidz(index) = (Psi(index) - Psi(index - 1)) / h;
                }
                else
                {
                    dPsidz(index) = (Psi(index + 1) - Psi(index - 1)) / (2 * h);
                }
                index++;
            }
        }
        // // 打印phi和dphidz
        // std::cout << "mode:" << mode << std::endl;
        // std::cout << "phi:\n" << eigen.phi.row(mode).real().transpose() << std::endl;

        // std::cout << "dphidz:\n" << eigen.dphidz.row(mode).real().transpose() << std::endl;

        for (int isz = 0; isz < params.Pos->NSz; isz++)
        {
            if (params.Pos->Sz(isz) > params.SSP[iprof].z(params.SSP[iprof].z.size() - 1))
            {
                // TODO处理海底
                eigen.PsiS(mode, isz) = 0.0;
                eigen.dPsidzS(mode, isz) = 0.0;
            }
            else
            {
                int index = ISzTab(isz);
                eigen.PsiS(mode, isz) = complex<double>(Psi(index)) + WTS(isz) * complex<double>(Psi(index + 1) - Psi(index));
                eigen.dPsidzS(mode, isz) = complex<double>(dPsidz(index)) + WTS(isz) * complex<double>(dPsidz(index + 1) - dPsidz(index));
            }
        }
        for (int irz = 0; irz < params.Pos->NRz; irz++)
        {
            if (params.Pos->Rz(irz) > params.SSP[iprof].z(params.SSP[iprof].z.size() - 1))
            {
                // TODO处理海底
                eigen.PsiR(mode, irz) = 0.0;
                eigen.dPsidzR(mode, irz) = 0.0;
            }
            else
            {
                int index = IRzTab(irz);
                eigen.PsiR(mode, irz) = complex<double>(Psi(index)) + WTR(irz) * complex<double>(Psi(index + 1) - Psi(index));
                eigen.dPsidzR(mode, irz) = complex<double>(dPsidz(index)) + WTR(irz) * complex<double>(dPsidz(index + 1) - dPsidz(index));
            }
        }
        for (int izz = 0; izz < NzTab; izz++)
        {
            if (zTab(izz) > params.SSP[iprof].z(params.SSP[iprof].z.size() - 1))
            {
                // TODO处理海底
                phiZ(mode, izz) = 0.0;
            }
            else
            {
                int index = IZzTab(izz);
                phiZ(mode, izz) = complex<double>(Psi(index)) + WTZ(izz) * complex<double>(Psi(index + 1) - Psi(index));
            }
        }

        // 输出phiZ到.mod中
        for (int izz = 0; izz < NzTab; izz++)
        {
            std::complex<float> phiZf = static_cast<std::complex<float>>(phiZ(mode, izz)); // 转换double为float
            // params.MODFile.write(reinterpret_cast<char *>(&phiZf), sizeof(std::complex<float>));
        }
    }
}

void Normalize(const size_t &iprof, const int &mode, const int &firstM, VectorXd &Phi, int &ITP, int &NTotal1, double &x, TridMtx &trid, const parameters &params, EigenParams &eigen)
{
    int iPower, modeCount = 0, j, j1, L, L1;
    double omega = (2 * pi * params.freqinfo->freq);
    double omega2 = SQ(omega);
    double x1, x2, DrhoDx = 0.0, DetaDx, SqNorm = 0.0, RN, ScaleFactor, Slow = 0.0; // rhoMedium, rho_omega_h2,
    double rhoMedium, rho_omega_h2;
    complex<double> Perturbation_k = 0.0, Del, fTop1, gTop1, fTop2, gTop2, fBot1, gBot1, fBot2, gBot2;
    bool isTop = false, isComplex = false;
    if (params.HSTop[iprof].BC == BC_Mode::MODE_A_Half_space)
    {
        complex<double> cp2 = SQ(params.HSTop[iprof].cp);
        Del = I1D * imag(sqrt(complex<double>(x - omega2 / cp2)));
        Perturbation_k -= Del * SQ(Phi(0)) / params.HSTop[iprof].rho;
        Slow += SQ(Phi(0)) / (2.0 * sqrt(x - real(omega2 / cp2))) / (params.HSTop[iprof].rho * real(cp2));
    }
    else if (params.HSTop[iprof].BC == BC_Mode::MODE_F_File || params.HSTop[iprof].BC == BC_Mode::MODE_P_Precomputed)
    {
        isTop = true;
        BCImpedance(iprof, x, isTop, fTop1, gTop1, iPower, isComplex, trid,
                    params, modeCount);
        isComplex = true;
        BCImpedance(iprof, x, isTop, fTop2, gTop2, iPower, isComplex, trid,
                    params, modeCount);
        Del = fTop2 / gTop2 - fTop1 / gTop1;
        Perturbation_k -= Del * SQ(Phi(0));
    }
    L = trid.Loc(params.SSP[iprof].FirstAcoustic) - 1;
    j = 0;
    for (int im = params.SSP[iprof].FirstAcoustic; im <= params.SSP[iprof].LastAcoustic; im++)
    {
        // Compute contribution from the volume
        L += 1;
        rhoMedium = trid.rho(L);
        rho_omega_h2 = rhoMedium * omega2 * SQ(trid.h(im));
        // top interface
        SqNorm += 0.5 * trid.h(im) * SQ(Phi(j)) / rhoMedium;
        Slow += 0.5 * trid.h(im) * (trid.B1(L) + 2.0) * SQ(Phi(j)) / rho_omega_h2;
        Perturbation_k += 0.5 * trid.h(im) * I1D * trid.B1C(L) * SQ(Phi(j)) / rhoMedium;

        // medium
        L1 = L + 1;
        L += trid.N(im) - 1;
        j1 = j + 1;
        j += trid.N(im) - 1;

        SqNorm += trid.h(im) * Phi.segment(j1, j - j1 + 1).array().square().sum() / rhoMedium;
        Slow += trid.h(im) * trid.B1.segment(j1, j - j1 + 1).array().square().sum() / rho_omega_h2;
        Perturbation_k += trid.h(im) * I1D * (trid.B1C.segment(L1, L - L1 + 1).array() * Phi.segment(j1, j - j1 + 1).array().square()).sum() / rhoMedium;

        // 底部边界：更新索引
        L += 1;
        j += 1;

        // 累加平方范数（SqNorm）
        SqNorm += 0.5 * trid.h(im) * (Phi(j) * Phi(j)) / rhoMedium;

        // 累加慢度项（Slow）
        Slow += 0.5 * trid.h(im) * (trid.B1(L) + 2.0) * SQ(Phi(j)) / rho_omega_h2;

        // 累加底部接口的扰动项（第二个Perturbation_k更新）
        Perturbation_k += 0.5 * trid.h(im) * I1D * trid.B1C(L) * SQ(Phi(j)) / rhoMedium;
    }

    if (params.HSBot[iprof].BC == BC_Mode::MODE_A_Half_space || params.HSBot[iprof].BC == BC_Mode::MODE_F_File || params.HSBot[iprof].BC == BC_Mode::MODE_P_Precomputed)
    {
        isTop = false;
        BCImpedance(iprof, x, isTop, fBot1, gBot1, iPower, isComplex, trid,
                    params, modeCount);
        isComplex = true;
        BCImpedance(iprof, x, isTop, fBot2, gBot2, iPower, isComplex, trid,
                    params, modeCount);
        Del = fBot2 / gBot2 - fBot1 / gBot1;
        Perturbation_k -= Del * SQ(Phi(j));
    }

    // Compute derivative of top admitance
    x1 = 0.9999999 * x;
    x2 = 1.0000001 * x;
    isTop = true;
    BCImpedance(iprof, x1, isTop, fTop1, gTop1, iPower, isComplex, trid,
                params, modeCount);
    BCImpedance(iprof, x2, isTop, fTop2, gTop2, iPower, isComplex, trid,
                params, modeCount);
    if (gTop1 != 0.0)
        DrhoDx = real((fTop2 / gTop2 - fTop1 / gTop1)) / (x2 - x1);

    // Compute derivative of bottom admitance
    isTop = false;
    BCImpedance(iprof, x1, isTop, fBot1, gBot1, iPower, isComplex, trid,
                params, modeCount);
    BCImpedance(iprof, x2, isTop, fBot2, gBot2, iPower, isComplex, trid,
                params, modeCount);
    if (gBot1 != 0.0)
        DetaDx = real((fBot2 / gBot2 - fBot1 / gBot1)) / (x2 - x1);

    // Scale the mode
    RN = SqNorm - DrhoDx * SQ(Phi(0)) + DetaDx * SQ(Phi(NTotal1 - 1));

    if (RN <= 0.0)
    {
        RN = -RN;
    }

    ScaleFactor = 1.0 / sqrt(RN);
    if (Phi(ITP) < 0.0)
        ScaleFactor = -ScaleFactor;

    Phi *= ScaleFactor;
    // std::cout << "Phi: " << Phi << std::endl;
    Slow = SQ(ScaleFactor) * Slow * omega / sqrt(x);
    Perturbation_k *= SQ(ScaleFactor);
    eigen.VG(mode) = 1 / Slow;
    ScatterLoss(iprof, mode, Perturbation_k, Phi, x, trid, params, eigen);
}

void ScatterLoss(const size_t &iprof, const int &mode, complex<double> &Perturbation_k, VectorXd &Psi, double &x, TridMtx &trid, const parameters &params, EigenParams &eigen)
{
    double omega = 2 * pi * params.freqinfo->freq;
    double omega2 = SQ(omega);
    double rho1, rho2, rhoInside, h2, sigma;
    complex<double> eta1Sq, eta2Sq, U, PhiC; // kx = sqrt(x),
    int j = 0, L = trid.Loc(params.SSP[iprof].FirstAcoustic);
    for (int im = params.SSP[iprof].FirstAcoustic - 1; im <= params.SSP[iprof].LastAcoustic; im++) // Loop over media
    {
        // Calculate rho1, eta1Sq, Phi, U
        if (im == params.SSP[iprof].FirstAcoustic - 1) // Top properties
        {
            switch (params.HSTop[iprof].BC)
            {
            case BC_Mode::MODE_A_Half_space:
                rho1 = params.HSTop[iprof].rho;
                eta1Sq = x - omega2 / SQ(params.HSTop[iprof].cp);
                U = sqrt(eta1Sq) * Psi(0) / rho1;
                break;
            case BC_Mode::MODE_V_Vacuum:
                rho1 = 1e-9;
                eta1Sq = 1.0;
                rhoInside = trid.rho(trid.Loc(params.SSP[iprof].FirstAcoustic));
                U = Psi(1) / trid.h(params.SSP[iprof].FirstAcoustic) / rhoInside;
                break;
            case BC_Mode::MODE_R_Rigid:
                rho1 = 1e9;
                eta1Sq = 1.0;
                U = 0.0;
                break;
            default:
                rho1 = 0.0;
                eta1Sq = 0.0;
                U = 0.0;
                break;
            }
            sigma = params.HSTop[iprof].sigma;
        }
        else
        {
            h2 = SQ(trid.h(im));
            j = j + trid.N(im);
            L = trid.Loc(im) + trid.N(im);
            rho1 = trid.rho(L);
            eta1Sq = (2.0 + trid.B1(L)) / h2 - x;
            U = (-Psi(j - 1) - 0.5 * (trid.B1(L) - h2 * x) * Psi(j)) / (trid.h(im) * rho1);
            sigma = params.SSP[iprof].sigma[im];
        }

        // Calculate rho2, eta2
        if (im == params.SSP[iprof].LastAcoustic) // Bottom properties
        {
            switch (params.HSBot[iprof].BC)
            {
            case BC_Mode::MODE_A_Half_space: // Acousto-elastic
                rho2 = params.HSBot[iprof].rho;
                eta2Sq = omega2 / params.HSBot[iprof].rho - x;
                break;
            case BC_Mode::MODE_V_Vacuum: // Vacuum
                rho2 = 1e-9;
                eta2Sq = 1.0;
                break;
            case BC_Mode::MODE_R_Rigid: // Rigid
                rho2 = 1e9;
                eta2Sq = 1.0;
                break;
            default: // Tabulated
                rho2 = 0.0;
                eta2Sq = 0.0;
                break;
            }
            sigma = params.HSBot[iprof].sigma;
        }
        else
        {
            rho2 = trid.rho(L + 1);
            eta2Sq = (2.0 + trid.B1(L + 1)) / SQ(trid.h(im + 1)) - x;
        }
        PhiC = Psi(j);
        complex<double> kup = KupIng(sigma, eta1Sq, rho1, eta2Sq, rho2, PhiC, U);
        // // 打印kup
        // cout << "kup: \n"
        //      << kup << endl;
        Perturbation_k += kup;
    }
    // 打印Perturbation_k
    // cout << "Perturbation_k: \n"
    //      << Perturbation_k << endl;
    eigen.k(mode) = Perturbation_k;
}