#include "solve.h"

// ERROUT函数声明（假设在其他文件中定义）
void ERROUT() {

};

void SolveEp(int &iset, const int &NSets, EigenParams &eigen, EigenFunction &eigenfun, KrakenMatrix &kramtrx, parameters &params, double &Error)
{
    double omega2 = SQ(2 * pi * params.freqinfo->freq);
    int iprof = 0;
    if (iprof > 0 && iset < 2 && params.modeType == ModeType::Couple)
    {
        Solve3();
    }
    else if ((iset < 2) && (params.NMedia <= params.LastAcoustic - params.FirstAcoustic + 1))
    {
        Solve1(iset, NSets, eigen, kramtrx, params);
    }
    else
    {
        Solve2(iset, eigen, kramtrx, params);
    }
    eigen.Extrap.resize(NSets * eigen.M);
    int start_idx = iset * eigen.M;
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
        VectorSolve(kramtrx, params, eigen, eigenfun, NTotal, NTotal1);
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
    int iPower = 0, NTotal, mode = 0; // NzTab = 0,
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

void Solve2(int &iset, EigenParams &eigen, KrakenMatrix &kramtrx, parameters &params)
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
                P(i) = eigen.EVMat(i * eigen.M + mode);
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
                x = P(1);
            }
        }
        Tolerance = abs(x) * kramtrx.B1.size() * pow(10.0, (1.0 - std::numeric_limits<double>::digits10));
        ZSecantX(x, Tolerance, Iteration, MaxIteration, iset, mode, Delta, iPower, kramtrx, params, eigen.EVMat, iscountm, modeCount, ErrorMessage, FUNCT);
        eigen.EVMat(iset * eigen.M + mode) = x;
        if (omega2 / SQ(params.Chigh) > x)
        {
            eigen.M--;
            return;
        }
    }
}

void Solve3(int &iset, const int &NSets, EigenParams &eigen, KrakenMatrix &kramtrx, parameters &params)
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

    FUNCT(iset, mode, xMin, Delta, iPower, kramtrx, params, eigen.EVMat, iscountm, modeCount);  
    int M = modeCount;

    for (int modeIdx = 0; modeIdx < M; ++modeIdx)
    {
        x = eigen.EVMat(iset * M + modeIdx);
        Tolerance = std::abs(x) * std::pow(10.0, 2.0 - std::numeric_limits<double>::digits10);
        ZSecantX(x, Tolerance, IT, MaxIT, iset, mode, Delta, iPower, kramtrx, params, eigen.EVMat, iscountm, modeCount, ErrorMessage, FUNCT);
        
        if (!ErrorMessage.empty())
        {
            // 输出警告信息
        }

        eigen.EVMat(iset * M + modeIdx) = x;

        if (omega2 / SQ(params.Chigh) > x)
        {
            eigen.M = modeIdx;  // 调整为当前索引
            return;
        }
    }
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

void VectorSolve(KrakenMatrix &kramtrx, parameters &params, EigenParams &eigen, EigenFunction &eigenfun,
                 int &NTotal, int &NTotal1)
{
    int j = 0, NzTab, iPower, modeCount = 0, L, ITP, IErr;
    double h_rho, x, xh2;
    VectorXd z(NTotal1), Phi(NTotal1), d(NTotal1), e(NTotal1 + 1), zTab, WTS, WTR;
    VectorXi ISzTab, IRzTab, Ix, Iy;
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
    // MergeVectors(params.Pos->Sz, params.Pos->Rz, zTab, NzTab, Ix, Iy);
    // 对声源深度和接收深度分别进行处理，不在需要合并

    eigenfun.phiS.resize(eigen.M, params.Pos->Sz.size());
    eigenfun.phiR.resize(eigen.M, params.Pos->Rz.size());

    WTS.resize(params.Pos->NSz);
    WTR.resize(params.Pos->NRz);
    ISzTab.resize(params.Pos->NSz);
    IRzTab.resize(params.Pos->NRz);

    Weight_dble(z, NTotal1, params.Pos->Sz, params.Pos->NSz, WTS, ISzTab);
    Weight_dble(z, NTotal1, params.Pos->Rz, params.Pos->NRz, WTR, IRzTab);

    // cout << "ISzTab: " << ISzTab.transpose() << endl;
    // cout << "IRzTab: " << IRzTab.transpose() << endl;
    // cout << "WTR:" << WTR.transpose() << endl; 

    // TODO 写入mod的头

    for (int mode = 0; mode < eigen.M; mode++)
    {
        x = eigen.EVMat(mode);
        BCImpedance(x, isTop, fTop, gTop, iPower, isComplex, kramtrx,
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
            xh2 = x * SQ(params.mesh.h(params.FirstAcoustic));
            h_rho = params.mesh.h(params.FirstAcoustic) * kramtrx.rho(L);
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
        BCImpedance(x, isTop, fBot, gBot, iPower, isComplex, kramtrx,
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
        }
        else
        {
            Normalize(mode, Phi, ITP, NTotal1, x, kramtrx, params, eigen);
        }

        for (int isz = 0; isz < params.Pos->NSz; isz++)
        {
            int index = ISzTab(isz);
            eigenfun.phiS(mode, isz) = complex<double>(Phi(index)) + WTS(isz) * complex<double>(Phi(index + 1) - Phi(index));
        }
        for (int irz = 0; irz < params.Pos->NRz; irz++)
        {
            int index = IRzTab(irz);
            eigenfun.phiR(mode, irz) = complex<double>(Phi(index)) + WTR(irz) * complex<double>(Phi(index + 1) - Phi(index));
        }
        // std::cout << "mode:" << mode << std::endl;
        // std::cout << "Phi: \n" << Phi << std::endl;
    }
}

void Normalize(int &mode, VectorXd &Phi, int &ITP, int &NTotal1, double &x, KrakenMatrix &kramtrx, parameters &params, EigenParams &eigen)
{
    int iPower, modeCount = 0, j, j1, L, L1;
    double omega = (2 * pi * params.freqinfo->freq);
    double omega2 = SQ(omega);
    double x1, x2, DrhoDx = 0.0, DetaDx, SqNorm = 0.0, RN, ScaleFactor, Slow = 0.0; // rhoMedium, rho_omega_h2,
    double rhoMedium, rho_omega_h2;
    complex<double> Perturbation_k = 0.0, Del, fTop1, gTop1, fTop2, gTop2, fBot1, gBot1, fBot2, gBot2;
    bool isTop = false, isComplex = false;
    if (params.HSTop.BC == BC_Mode::MODE_A_Half_space)
    {
        complex<double> cp2 = SQ(params.HSTop.cp);
        Del = I1D * imag(sqrt(complex<double>(x - omega2 / cp2)));
        Perturbation_k -= Del * SQ(Phi(0)) / params.HSTop.rho;
        Slow += SQ(Phi(0)) / (2.0 * sqrt(x - real(omega2 / cp2))) / (params.HSTop.rho * real(cp2));
    }
    else if (params.HSTop.BC == BC_Mode::MODE_F_File || params.HSTop.BC == BC_Mode::MODE_P_Precomputed)
    {
        isTop = true;
        BCImpedance(x, isTop, fTop1, gTop1, iPower, isComplex, kramtrx,
                    params, modeCount);
        isComplex = true;
        BCImpedance(x, isTop, fTop2, gTop2, iPower, isComplex, kramtrx,
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

    if (params.HSBot.BC == BC_Mode::MODE_A_Half_space || params.HSBot.BC == BC_Mode::MODE_F_File || params.HSBot.BC == BC_Mode::MODE_P_Precomputed)
    {
        isTop = false;
        BCImpedance(x, isTop, fBot1, gBot1, iPower, isComplex, kramtrx,
                    params, modeCount);
        isComplex = true;
        BCImpedance(x, isTop, fBot2, gBot2, iPower, isComplex, kramtrx,
                    params, modeCount);
        Del = fBot2 / gBot2 - fBot1 / gBot1;
        Perturbation_k -= Del * SQ(Phi(j));
    }

    // Compute derivative of top admitance
    x1 = 0.9999999 * x;
    x2 = 1.0000001 * x;
    isTop = true;
    BCImpedance(x1, isTop, fTop1, gTop1, iPower, isComplex, kramtrx,
                params, modeCount);
    BCImpedance(x2, isTop, fTop2, gTop2, iPower, isComplex, kramtrx,
                params, modeCount);
    if (gTop1 != 0.0)
        DrhoDx = real((fTop2 / gTop2 - fTop1 / gTop1)) / (x2 - x1);

    // Compute derivative of bottom admitance
    isTop = false;
    BCImpedance(x1, isTop, fBot1, gBot1, iPower, isComplex, kramtrx,
                params, modeCount);
    BCImpedance(x2, isTop, fBot2, gBot2, iPower, isComplex, kramtrx,
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
    Slow = SQ(ScaleFactor) * Slow * omega / sqrt(x);
    Perturbation_k *= SQ(ScaleFactor);
    eigen.VG(mode) = 1 / Slow;
    ScatterLoss(mode, Perturbation_k, Phi, x, kramtrx, params, eigen);
}

void ScatterLoss(int &mode, complex<double> &Perturbation_k, VectorXd &Phi, double &x, KrakenMatrix &kramtrx, parameters &params, EigenParams &eigen)
{
    double omega = 2 * pi * params.freqinfo->freq;
    double omega2 = SQ(omega);
    double rho1, rho2, rhoInside, h2;
    complex<double> eta1Sq, eta2Sq, U, PhiC; // kx = sqrt(x),
    int j = 0, L = params.mesh.Loc(params.FirstAcoustic);
    for (int im = params.FirstAcoustic - 1; im <= params.LastAcoustic; im++) // Loop over media
    {
        // Calculate rho1, eta1Sq, Phi, U
        if (im == params.FirstAcoustic - 1) // Top properties
        {
            switch (params.HSTop.BC)
            {
            case BC_Mode::MODE_A_Half_space:
                rho1 = params.HSTop.rho;
                eta1Sq = x - omega2 / SQ(params.HSTop.cp);
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
        }
        else
        {
            h2 = SQ(params.mesh.h(im));
            j = j + params.mesh.N(im);
            L = params.mesh.Loc(im) + params.mesh.N(im);
            rho1 = kramtrx.rho(L);
            eta1Sq = (2.0 + kramtrx.B1(L)) / h2 - x;
            U = (-Phi(j - 1) - 0.5 * (kramtrx.B1(L) - h2 * x) * Phi(j)) / (params.mesh.h(im) * rho1);
        }

        // Calculate rho2, eta2
        if (im == params.LastAcoustic) // Bottom properties
        {
            switch (params.HSBot.BC)
            {
            case BC_Mode::MODE_A_Half_space: // Acousto-elastic
                rho2 = params.HSBot.rho;
                eta2Sq = omega2 / params.HSBot.rho - x;
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
        }
        else
        {
            rho2 = kramtrx.rho(L + 1);
            eta2Sq = (2.0 + kramtrx.B1(L + 1)) / SQ(params.mesh.h(im + 1)) - x;
        }
        PhiC = Phi(j);
        complex<double> kup = KupIng(params.SSP[im + 1].sigma, eta1Sq, rho1, eta2Sq, rho2, PhiC, U);
        Perturbation_k += kup;
    }
    eigen.k(mode) = Perturbation_k;
}