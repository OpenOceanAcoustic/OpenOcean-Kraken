#include "solve.h"

// ERROUT函数声明（假设在其他文件中定义）
void ERROUT() {

};

void SolveEp(int &iset, size_t iprof, const int &NSets, EigenParams &eigen, EigenFunction &eigenfun, KrakenMatrix &kramtrx, parameters &params, double &Error)
{
    double omega2 = SQ(2 * pi * params.freqinfo->freq);
    if (iprof > 0 && iset < 2 && params.modeType == ModeType::Couple)
    {
        Solve3(iset, iprof, eigen, kramtrx, params);
    }
    else if ((iset < 2) && (params.SSP[iprof].NMedia <= params.LastAcoustic - params.FirstAcoustic + 1))
    {
        Solve1(iset, iprof, NSets, eigen, kramtrx, params);
    }
    else
    {
        Solve2(iset, iprof, eigen, kramtrx, params);
    }
    if (iset == 0)
    {
        eigen.firstM = eigen.M;
        eigen.Extrap.resize(NSets * eigen.M);
    }

    int start_idx = iset * eigen.firstM;
    eigen.Extrap.segment(start_idx, eigen.M) = eigen.EVMat.segment(start_idx, eigen.M);

    // std::cout << "iset: " << iset << " \n"
    //           << eigen.EVMat.segment(start_idx, eigen.M) << std::endl;

    // 查找满足条件的最小位置
    // Extrap(1, 1:M)对应索引为0到M-1
    VectorXd Ex1 = eigen.Extrap.segment(0, eigen.M);
    double threshold = omega2 / SQ(params.Chigh);

    int Min_Loc = 0;

    while (Min_Loc < eigen.M && Ex1(Min_Loc) > threshold)
    {
        Min_Loc++;
    }
    if (Min_Loc != 0)
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
        VectorSolve(iprof, kramtrx, params, eigen, eigenfun, NTotal, NTotal1);
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
void Solve1(int &iset, size_t iprof, const int &NSets, EigenParams &eigen, KrakenMatrix &kramtrx, parameters &params)
{
    int iPower = 0, NTotal, mode = 0; // NzTab = 0,
    double x, x1, x2, xMin, xMax, Eps, Delta;
    bool iscountm = true;
    std::string ErrorMessage;
    VectorXd xL, xR;
    double omega2 = SQ(2 * pi * params.freqinfo->freq);

    int modeCount = 0;

    // 确定模态数量
    xMin = 1.00001 * omega2 / SQ(params.Chigh);

    FUNCT(iset, iprof, mode, xMin, Delta, iPower, kramtrx, params, eigen.EVMat, eigen.firstM, iscountm, modeCount);
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
    FUNCT(iset, iprof, mode, xMax, Delta, iPower, kramtrx, params, eigen.EVMat, eigen.firstM, iscountm, modeCount);

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
    for (int i = params.FirstAcoustic; i <= params.LastAcoustic; ++i)
    {
        NTotal += params.mesh.N(i);
    }

    if (M > NTotal / 5)
    {
        // 输出警告信息
        std::cout << "Approximate number of modes = " << M << " Warning in KRAKEN - Solve1 : Mesh too coarse to sample the modes adequately" << std::endl;
    }

    Bisection(iset, iprof, mode, xMin, xMax, xL, xR, kramtrx, params, eigen); // 初始化上下边界

    // 使用ZBRENT精化每个本征值
    iscountm = false;

    for (int modeIdx = 0; modeIdx < M; ++modeIdx)
    {
        x1 = xL(modeIdx);
        x2 = xR(modeIdx);
        Eps = std::abs(x2) * std::pow(10.0, 2.0 - std::numeric_limits<double>::digits10);
        ZBRENTX(x, x1, x2, Eps, iset, iprof, mode, Delta, iPower, kramtrx, params, eigen.EVMat, eigen.firstM, iscountm, modeCount,
                ErrorMessage, FUNCT); // Brent求根法

        if (!ErrorMessage.empty())
        {
            // 输出警告信息
        }

        eigen.EVMat(iset * eigen.firstM + modeIdx) = x;
    }
}

void Solve2(int &iset, size_t iprof, EigenParams &eigen, KrakenMatrix &kramtrx, parameters &params)
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
                        x1 = SQ(params.mesh.hV(j));
                        x2 = SQ(params.mesh.hV(j + ii + 1));
                        P(j) = ((SQ(params.mesh.hV(iset)) - x2) * P(j) -
                                (SQ(params.mesh.hV(iset)) - x1) * P(j + 1)) /
                               (x1 - x2);
                    }
                }
                x = P(0);
            }
        }
        Tolerance = abs(x) * kramtrx.B1.size() * pow(10.0, (1.0 - std::numeric_limits<double>::digits10));
        ZSecantX(x, Tolerance, Iteration, MaxIteration, iset, iprof, mode, Delta, iPower, kramtrx, params, eigen.EVMat, eigen.firstM, iscountm, modeCount, ErrorMessage, FUNCT);
        eigen.EVMat(iset * eigen.firstM + mode) = x;
        if (omega2 / SQ(params.Chigh) > x)
        {
            eigen.M = mode;
            return;
        }
    }
}

void Solve3(int &iset, size_t iprof, EigenParams &eigen, KrakenMatrix &kramtrx, parameters &params)
{
    int IT, MaxIT, iPower = 0, mode = 0; // NzTab = 0,
    double x, xMin, Tolerance, Delta;
    std::string ErrorMessage;

    bool iscountm = false;

    MaxIT = 500;

    double omega2 = SQ(2 * pi * params.freqinfo->freq);

    int modeCount = 0;

    // 确定模态数量
    xMin = 1.00001 * omega2 / SQ(params.Chigh);

    FUNCT(iset, iprof, mode, xMin, Delta, iPower, kramtrx, params, eigen.EVMat, eigen.firstM, iscountm, modeCount);
    int M = modeCount;

    for (int modeIdx = 0; modeIdx < M; ++modeIdx)
    {
        x = eigen.EVMat(iset * eigen.firstM + modeIdx);
        Tolerance = std::abs(x) * std::pow(10.0, 2.0 - std::numeric_limits<double>::digits10);
        ZSecantX(x, Tolerance, IT, MaxIT, iset, iprof, mode, Delta, iPower, kramtrx, params, eigen.EVMat, eigen.firstM, iscountm, modeCount, ErrorMessage, FUNCT);

        if (!ErrorMessage.empty())
        {
            // 输出警告信息
        }

        eigen.EVMat(iset * eigen.firstM + modeIdx) = x;

        if (omega2 / SQ(params.Chigh) > x)
        {
            eigen.M = modeIdx; // 调整为当前索引
            return;
        }
    }
}

// FUNCT函数：计算色散关系
void FUNCT(int &iset, size_t iprof, int &mode, double &x, double &Delta, int &iPower, KrakenMatrix &kramtrx,
           parameters &params, VectorXd &EVMat, const int &firstM, bool coutmodes, int &modeCount)
{
    int iPowerBot;
    double f, g;
    std::complex<double> fTop, gTop, fBot, gBot;
    bool isTop = false, isComplex = false;

    modeCount = 0;

    // 调用BCImpedance计算底部阻抗
    BCImpedance(iprof, x, isTop, fTop, gTop, iPower, isComplex, kramtrx,
                params, modeCount);

    f = std::real(fTop);
    g = std::real(gTop);

    // 穿过声学层
    AcousticLayers(x, f, g, iPower, kramtrx, params, coutmodes, modeCount);

    isTop = true;
    // 调用BCImpedance计算顶部阻抗
    BCImpedance(iprof, x, isTop, fBot, gBot, iPowerBot, isComplex, kramtrx,
                params, modeCount);

    Delta = std::real(f * std::real(gBot) - g * std::real(fBot));
    iPower = iPower + iPowerBot;

    if (g * Delta > 0.0)
    {
        modeCount++;
    }

    // 减去之前的根
    if (mode > 0 && params.SSP[iprof].NMedia > params.LastAcoustic - params.FirstAcoustic + 1)
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
void Bisection(int &iset, size_t iprof, int &mode, double xMin, double xMax, VectorXd &xL, VectorXd &xR,
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
    FUNCT(iset, iprof, mode, xMax, Delta, iPower, kramtrx, params, eigen.EVMat, eigen.firstM, true, NZer1);

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
                FUNCT(iset, iprof, modeIdx, x, Delta, iPower, kramtrx, params, eigen.EVMat, eigen.firstM, true, modeCount);
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

void VectorSolve(size_t iprof, KrakenMatrix &kramtrx, parameters &params, EigenParams &eigen, EigenFunction &eigenfun,
                 int &NTotal, int &NTotal1)
{
    int j = 0, NzTab, iPower, modeCount = 0, L, ITP, IErr;
    double h_rho, x, xh2;
    VectorXd z(NTotal1), Phi(NTotal1), d(NTotal1), e(NTotal1 + 1), zTab, WTS, WTR, WTZ;
    VectorXi ISzTab, IRzTab, IZzTab, Ix, Iy;
    VectorXcd PhiTab;
    complex<double> fTop, gTop, fBot, gBot;
    bool isTop = true, isComplex = false;
    z(0) = params.SSP->z(params.FirstAcoustic);
    for (int im = params.FirstAcoustic; im <= params.LastAcoustic; im++)
    {
        h_rho = params.mesh.h(im) * kramtrx.rho(params.mesh.Loc(im));
        for (int ii = 1; ii <= params.mesh.N(im); ii++)
        {
            e(j + ii) = 1.0 / h_rho;
            z(j + ii) = z(j) + params.mesh.h(im) * ii;
        }
        j = j + params.mesh.N(im);
    }
    e(NTotal1) = 1.0 / h_rho;
    // 为了方便输出.mod，这里还是进行合并
    MergeVectors(params.Pos->Sz, params.Pos->Rz, zTab, NzTab, Ix, Iy);
    MatrixXcd phiZ(eigen.firstM, NzTab);

    // 对声源深度和接收深度分别进行处理，不在需要合并

    eigenfun.phi.resize(eigen.firstM, NTotal1);
    eigenfun.phiS.resize(eigen.firstM, params.Pos->Sz.size());
    eigenfun.phiR.resize(eigen.firstM, params.Pos->Rz.size());
    eigenfun.dphidz.resize(eigen.firstM, NTotal1);
    eigenfun.dphidzS.resize(eigen.firstM, params.Pos->Sz.size());
    eigenfun.dphidzR.resize(eigen.firstM, params.Pos->Rz.size());

    WTS.resize(params.Pos->NSz);
    WTR.resize(params.Pos->NRz);
    WTZ.resize(NzTab);

    ISzTab.resize(params.Pos->NSz);
    IRzTab.resize(params.Pos->NRz);
    IZzTab.resize(NzTab);

    Weight_dble(z, NTotal1, params.Pos->Sz, params.Pos->NSz, WTS, ISzTab);
    Weight_dble(z, NTotal1, params.Pos->Rz, params.Pos->NRz, WTR, IRzTab);
    Weight_dble(z, NTotal1, zTab, NzTab, WTZ, IZzTab);

    // cout << "ISzTab: " << ISzTab.transpose() << endl;
    // cout << "IRzTab: " << IRzTab.transpose() << endl;
    // cout << "WTR:" << WTR.transpose() << endl;

    // TODO 写入mod的头
    string filename = "test";
    params.MODFile.open(filename + ".mod", std::ios::binary);
    if (!params.MODFile.is_open())
    {
        std::cerr << "无法打开MOD文件: " << filename << ".mod" << std::endl;
        return;
    }

    int ifreq = 0;
    if (ifreq == 0 && iprof == 0)
    {
        eigen.LRecordLength = std::max(2 * params.freqinfo->Nfreq, std::max(2 * NzTab, std::max(32, 3 * (params.LastAcoustic - params.FirstAcoustic + 1))));
    }
    if (ifreq == 0)
    {
        eigen.IRecProfile = 0;
        params.MODFile.seekp(eigen.IRecProfile * 4 * eigen.LRecordLength, std::ios::beg);
        params.MODFile.write(reinterpret_cast<char *>(&eigen.LRecordLength), sizeof(int));
        // 固定title为80个char
        char title[80] = {0};
        strncpy(title, params.Title.c_str(), params.Title.size());
        params.MODFile.write(title, sizeof(title));
        params.MODFile.write(reinterpret_cast<char *>(&params.freqinfo->Nfreq), sizeof(int));
        int numAcoustic = params.LastAcoustic - params.FirstAcoustic + 1;
        params.MODFile.write(reinterpret_cast<char *>(&numAcoustic), sizeof(int));
        params.MODFile.write(reinterpret_cast<char *>(&NzTab), sizeof(int));
        params.MODFile.write(reinterpret_cast<char *>(&NzTab), sizeof(int));

        params.MODFile.seekp((eigen.IRecProfile + 1) * 4 * eigen.LRecordLength, std::ios::beg);
        for (int im = params.FirstAcoustic; im <= params.LastAcoustic; im++)
        {
            params.MODFile.write(reinterpret_cast<char *>(&params.mesh.N(im)), sizeof(int));
            params.MODFile.write(reinterpret_cast<char *>(&params.SSP->Material), sizeof(params.SSP->Material)); // fortran是字符串，这里是整数。而且只有一层，应该是params.SSP->Material[im]?
        }
        params.MODFile.seekp((eigen.IRecProfile + 2) * 4 * eigen.LRecordLength, std::ios::beg);
        for (int im = params.FirstAcoustic; im <= params.LastAcoustic; im++)
        {
            float depth = static_cast<float>(params.SSP[iprof].depth[im]); // 转换double为float
            params.MODFile.write(reinterpret_cast<char *>(&depth), sizeof(float));
            float rho = kramtrx.rho(params.mesh.Loc(im)); // 转换double为float
            params.MODFile.write(reinterpret_cast<char *>(&rho), sizeof(float));
        }
        params.MODFile.seekp((eigen.IRecProfile + 3) * 4 * eigen.LRecordLength, std::ios::beg);
        params.MODFile.write(reinterpret_cast<char *>(params.freqinfo->freqvec.data()), params.freqinfo->Nfreq * sizeof(double));
        params.MODFile.seekp((eigen.IRecProfile + 4) * 4 * eigen.LRecordLength, std::ios::beg);
        // zTab转为float输出
        VectorXf zTabf = zTab.cast<float>();
        params.MODFile.write(reinterpret_cast<char *>(zTabf.data()), NzTab * sizeof(float));
        // for (int isz = 0; isz < params.Pos->NSz; isz++)
        // {
        //     params.MODFile.write(reinterpret_cast<char*>(&params.Pos->Sz(isz)), sizeof(double));
        // }
        // for (int irz = 0; irz < params.Pos->NRz; irz++)
        // {
        //     params.MODFile.write(reinterpret_cast<char*>(&params.Pos->Rz(irz)), sizeof(double));
        // }
        eigen.IRecProfile += 5;
    }
    params.MODFile.seekp((eigen.IRecProfile + 1) * 4 * eigen.LRecordLength, std::ios::beg);
    params.MODFile.write(reinterpret_cast<char *>(&params.HSTop[iprof].BC), sizeof(params.HSTop[iprof].BC)); // 这也是整数而不是字符串
    params.MODFile.write(reinterpret_cast<char *>(&params.HSTop[iprof].cp), sizeof(std::complex<double>));
    params.MODFile.write(reinterpret_cast<char *>(&params.HSTop[iprof].cs), sizeof(std::complex<double>));
    float HSToprho = static_cast<float>(params.HSTop[iprof].rho); // 转换double为float
    params.MODFile.write(reinterpret_cast<char *>(&HSToprho), sizeof(float));
    float SSPz0 = static_cast<float>(params.SSP[iprof].z(0)); // 转换double为float
    params.MODFile.write(reinterpret_cast<char *>(&SSPz0), sizeof(float));
    params.MODFile.write(reinterpret_cast<char *>(&params.HSBot[iprof].BC), sizeof(params.HSBot[iprof].BC)); // 这也是整数而不是字符串
    params.MODFile.write(reinterpret_cast<char *>(&params.HSBot[iprof].cp), sizeof(std::complex<double>));
    params.MODFile.write(reinterpret_cast<char *>(&params.HSBot[iprof].cs), sizeof(std::complex<double>));
    float HSBotrho = static_cast<float>(params.HSBot[iprof].rho); // 转换double为float
    params.MODFile.write(reinterpret_cast<char *>(&HSBotrho), sizeof(float));
    float SSPzn = static_cast<float>(params.SSP[iprof].z(params.SSP[iprof].z.size() - 1)); // 转换double为float
    params.MODFile.write(reinterpret_cast<char *>(&SSPzn), sizeof(float));

    for (int mode = 0; mode < eigen.firstM; mode++)
    {
        x = eigen.EVMat(mode);
        BCImpedance(iprof, x, isTop, fTop, gTop, iPower, isComplex, kramtrx,
                    params, modeCount);
        if (gTop == 0.0)
        {
            d(0) = 1.0;
            e(1) = 0.0;
        }
        else
        {
            L = params.mesh.Loc(params.FirstAcoustic);
            xh2 = x * SQ(params.mesh.h(params.FirstAcoustic));
            h_rho = params.mesh.h(params.FirstAcoustic) * kramtrx.rho(L);
            d(0) = (kramtrx.B1(L) - xh2) / h_rho / 2.0 + real(fTop / gTop);
        }

        ITP = NTotal;
        j = 0;
        L = params.mesh.Loc(params.FirstAcoustic);
        for (int im = params.FirstAcoustic; im <= params.LastAcoustic; im++)
        {
            xh2 = x * SQ(params.mesh.h(im));
            h_rho = params.mesh.h(im) * kramtrx.rho(L + 1);
            if (im >= params.FirstAcoustic + 1)
            {
                L++;
                d(j) = (d(j) + (kramtrx.B1(L) - xh2) / h_rho) / 2.0;
            }
            for (int ii = 0; ii < params.mesh.N(im); ii++)
            {
                j++;
                L++;
                d(j) = (kramtrx.B1(L) - xh2) / h_rho;
                if (kramtrx.B1(L) - xh2 + 2.0 > 0.0)
                    ITP = min(j, ITP);
            }
        }
        isTop = false;
        BCImpedance(iprof, x, isTop, fBot, gBot, iPower, isComplex, kramtrx,
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
        InverseIterationD(NTotal1, d, e, IErr, Phi);

        if (IErr != 0)
        {
            // 未收敛
            std::cout << "mode = " << mode << std::endl;
            std::cout << "Warning in KRAKEN - InverseIteration: Inverse iteration failed to converge" << std::endl;
            Phi = VectorXd::Zero(NTotal1); // zero out the errant eigenvector
        }
        else
        {
            Normalize(iprof, mode, eigen.firstM, Phi, ITP, NTotal1, x, kramtrx, params, eigen);
        }
        int index = 0;
        for (int im = params.FirstAcoustic; im <= params.LastAcoustic; ++im)
        {
            double h = params.mesh.h(im);
            int Nn = params.mesh.N(im);
            if (im == params.LastAcoustic)
            {
                Nn++;
            }
            for (int ii = 0; ii < Nn; ii++)
            {
                eigenfun.phi(mode, index) = Phi(index);
                if (index == 0)
                {
                    eigenfun.dphidz(mode, index) = (Phi(index + 1) - Phi(index)) / h;
                }
                else if (index == NTotal)
                {
                    eigenfun.dphidz(mode, index) = (Phi(index) - Phi(index - 1)) / h;
                }
                else
                {
                    eigenfun.dphidz(mode, index) = (Phi(index + 1) - Phi(index - 1)) / (2 * h);
                }
                index++;
            }
        }
        // // 打印phi和dphidz
        // std::cout << "mode:" << mode << std::endl;
        // std::cout << "phi:\n" << eigenfun.phi.row(mode).real().transpose() << std::endl;

        // std::cout << "dphidz:\n" << eigenfun.dphidz.row(mode).real().transpose() << std::endl;

        for (int isz = 0; isz < params.Pos->NSz; isz++)
        {
            if (params.Pos->Sz(isz) > params.SSP[iprof].z(params.SSP[iprof].z.size() - 1))
            {
                // TODO处理海底
                eigenfun.phiS(mode, isz) = 0.0;
                eigenfun.dphidzS(mode, isz) = 0.0;
            }
            else
            {
                int index = ISzTab(isz);
                eigenfun.phiS(mode, isz) = complex<double>(Phi(index)) + WTS(isz) * complex<double>(Phi(index + 1) - Phi(index));
                eigenfun.dphidzS(mode, isz) = complex<double>(eigenfun.dphidz(mode, index)) + WTS(isz) * complex<double>(eigenfun.dphidz(mode, index + 1) - eigenfun.dphidz(mode, index));
            }
        }
        for (int irz = 0; irz < params.Pos->NRz; irz++)
        {
            if (params.Pos->Rz(irz) > params.SSP[iprof].z(params.SSP[iprof].z.size() - 1))
            {
                // TODO处理海底
                eigenfun.phiR(mode, irz) = 0.0;
                eigenfun.dphidzR(mode, irz) = 0.0;
            }
            else
            {
                int index = IRzTab(irz);
                eigenfun.phiR(mode, irz) = complex<double>(Phi(index)) + WTR(irz) * complex<double>(Phi(index + 1) - Phi(index));
                eigenfun.dphidzR(mode, irz) = complex<double>(eigenfun.dphidz(mode, index)) + WTR(irz) * complex<double>(eigenfun.dphidz(mode, index + 1) - eigenfun.dphidz(mode, index));
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
                phiZ(mode, izz) = complex<double>(Phi(index)) + WTZ(izz) * complex<double>(Phi(index + 1) - Phi(index));
            }
        }
        // std::cout << "mode:" << mode << std::endl;
        // std::cout << "Phi: \n" << Phi << std::endl;
        // 写入Phi
        params.MODFile.seekp((eigen.IRecProfile + 2 + mode) * 4 * eigen.LRecordLength, std::ios::beg);

        // for (int isz = 0; isz < params.Pos->NSz; isz++)
        // {
        //     std::complex<float> phiS = static_cast<std::complex<float>>(eigenfun.phiS(mode, isz));  // 转换double为float
        //     MODFile.write(reinterpret_cast<char*>(&phiS), sizeof(std::complex<float>));
        // }
        // for (int irz = 0; irz < params.Pos->NRz; irz++)
        // {
        //     std::complex<float> phiR = static_cast<std::complex<float>>(eigenfun.phiR(mode, irz));  // 转换double为float
        //     MODFile.write(reinterpret_cast<char*>(&phiR), sizeof(std::complex<float>));
        // }
        // 输出phiZ到.mod中
        for (int izz = 0; izz < NzTab; izz++)
        {
            std::complex<float> phiZf = static_cast<std::complex<float>>(phiZ(mode, izz)); // 转换double为float
            params.MODFile.write(reinterpret_cast<char *>(&phiZf), sizeof(std::complex<float>));
        }
    }
}

void Normalize(size_t iprof, int &mode, int &firstM, VectorXd &Phi, int &ITP, int &NTotal1, double &x, KrakenMatrix &kramtrx, parameters &params, EigenParams &eigen)
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
        BCImpedance(iprof, x, isTop, fTop1, gTop1, iPower, isComplex, kramtrx,
                    params, modeCount);
        isComplex = true;
        BCImpedance(iprof, x, isTop, fTop2, gTop2, iPower, isComplex, kramtrx,
                    params, modeCount);
        Del = fTop2 / gTop2 - fTop1 / gTop1;
        Perturbation_k -= Del * SQ(Phi(0));
    }
    L = params.mesh.Loc(params.FirstAcoustic) - 1;
    j = 0;
    for (int im = params.FirstAcoustic; im <= params.LastAcoustic; im++)
    {
        // Compute contribution from the volume
        L += 1;
        rhoMedium = kramtrx.rho(L);
        rho_omega_h2 = rhoMedium * omega2 * SQ(params.mesh.h(im));
        // top interface
        SqNorm += 0.5 * params.mesh.h(im) * SQ(Phi(j)) / rhoMedium;
        Slow += 0.5 * params.mesh.h(im) * (kramtrx.B1(L) + 2.0) * SQ(Phi(j)) / rho_omega_h2;
        Perturbation_k += 0.5 * params.mesh.h(im) * I1D * kramtrx.B1C(L) * SQ(Phi(j)) / rhoMedium;

        // medium
        L1 = L + 1;
        L += params.mesh.N(im) - 1;
        j1 = j + 1;
        j += params.mesh.N(im) - 1;

        SqNorm += params.mesh.h(im) * Phi.segment(j1, j - j1 + 1).array().square().sum() / rhoMedium;
        Slow += params.mesh.h(im) * kramtrx.B1.segment(j1, j - j1 + 1).array().square().sum() / rho_omega_h2;
        Perturbation_k += params.mesh.h(im) * I1D * (kramtrx.B1C.segment(L1, L - L1 + 1).array() * Phi.segment(j1, j - j1 + 1).array().square()).sum() / rhoMedium;

        // 底部边界：更新索引
        L += 1;
        j += 1;

        // 累加平方范数（SqNorm）
        SqNorm += 0.5 * params.mesh.h(im) * (Phi(j) * Phi(j)) / rhoMedium;

        // 累加慢度项（Slow）
        Slow += 0.5 * params.mesh.h(im) * (kramtrx.B1(L) + 2.0) * SQ(Phi(j)) / rho_omega_h2;

        // 累加底部接口的扰动项（第二个Perturbation_k更新）
        Perturbation_k += 0.5 * params.mesh.h(im) * I1D * kramtrx.B1C(L) * SQ(Phi(j)) / rhoMedium;
    }

    if (params.HSBot[iprof].BC == BC_Mode::MODE_A_Half_space || params.HSBot[iprof].BC == BC_Mode::MODE_F_File || params.HSBot[iprof].BC == BC_Mode::MODE_P_Precomputed)
    {
        isTop = false;
        BCImpedance(iprof, x, isTop, fBot1, gBot1, iPower, isComplex, kramtrx,
                    params, modeCount);
        isComplex = true;
        BCImpedance(iprof, x, isTop, fBot2, gBot2, iPower, isComplex, kramtrx,
                    params, modeCount);
        Del = fBot2 / gBot2 - fBot1 / gBot1;
        Perturbation_k -= Del * SQ(Phi(j));
    }

    // Compute derivative of top admitance
    x1 = 0.9999999 * x;
    x2 = 1.0000001 * x;
    isTop = true;
    BCImpedance(iprof, x1, isTop, fTop1, gTop1, iPower, isComplex, kramtrx,
                params, modeCount);
    BCImpedance(iprof, x2, isTop, fTop2, gTop2, iPower, isComplex, kramtrx,
                params, modeCount);
    if (gTop1 != 0.0)
        DrhoDx = real((fTop2 / gTop2 - fTop1 / gTop1)) / (x2 - x1);

    // Compute derivative of bottom admitance
    isTop = false;
    BCImpedance(iprof, x1, isTop, fBot1, gBot1, iPower, isComplex, kramtrx,
                params, modeCount);
    BCImpedance(iprof, x2, isTop, fBot2, gBot2, iPower, isComplex, kramtrx,
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
    ScatterLoss(iprof, mode, Perturbation_k, Phi, x, kramtrx, params, eigen);
}

void ScatterLoss(size_t iprof, int &mode, complex<double> &Perturbation_k, VectorXd &Phi, double &x, KrakenMatrix &kramtrx, parameters &params, EigenParams &eigen)
{
    double omega = 2 * pi * params.freqinfo->freq;
    double omega2 = SQ(omega);
    double rho1, rho2, rhoInside, h2, sigma;
    complex<double> eta1Sq, eta2Sq, U, PhiC; // kx = sqrt(x),
    int j = 0, L = params.mesh.Loc(params.FirstAcoustic);
    for (int im = params.FirstAcoustic - 1; im <= params.LastAcoustic; im++) // Loop over media
    {
        // Calculate rho1, eta1Sq, Phi, U
        if (im == params.FirstAcoustic - 1) // Top properties
        {
            switch (params.HSTop[iprof].BC)
            {
            case BC_Mode::MODE_A_Half_space:
                rho1 = params.HSTop[iprof].rho;
                eta1Sq = x - omega2 / SQ(params.HSTop[iprof].cp);
                U = sqrt(eta1Sq) * Phi(0) / rho1;
                break;
            case BC_Mode::MODE_V_Vacuum:
                rho1 = 1e-9;
                eta1Sq = 1.0;
                rhoInside = kramtrx.rho(params.mesh.Loc(params.FirstAcoustic));
                U = Phi(1) / params.mesh.h(params.FirstAcoustic) / rhoInside;
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
            h2 = SQ(params.mesh.h(im));
            j = j + params.mesh.N(im);
            L = params.mesh.Loc(im) + params.mesh.N(im);
            rho1 = kramtrx.rho(L);
            eta1Sq = (2.0 + kramtrx.B1(L)) / h2 - x;
            U = (-Phi(j - 1) - 0.5 * (kramtrx.B1(L) - h2 * x) * Phi(j)) / (params.mesh.h(im) * rho1);
            sigma = params.SSP[iprof].sigma[im];
        }

        // Calculate rho2, eta2
        if (im == params.LastAcoustic) // Bottom properties
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
            rho2 = kramtrx.rho(L + 1);
            eta2Sq = (2.0 + kramtrx.B1(L + 1)) / SQ(params.mesh.h(im + 1)) - x;
        }
        PhiC = Phi(j);
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