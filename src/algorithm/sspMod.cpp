#include "sspMod.h"

void EvaluateSSP(SSPStructure &SSP)
{
    // double cxx, cyy, cxy, cxz, cyz;

    switch (SSP.Type)
    {
    case SSP_Mode::MODE_N_n2Linear:
        n2Linear(SSP);
        break;
    case SSP_Mode::MODE_C_cLinear:
        cLinear(SSP);
        break;
    case SSP_Mode::MODE_P_cPCHIP:
        cPCHIP(SSP);
        break;
    case SSP_Mode::MODE_S_cCubic:
        cCubic(SSP);
        break;
    // case SSP_Mode::MODE_A_Analytic:
    //     Analytic(cp, cs, rho_k, Medium, N1); // 需要实现Analytic函数
    //     break;
    default:
        // 报错
        std::cerr << "Unknown SSP type in EvaluateSSP" << std::endl;
        break;
    }
}

/**
 * @brief 根据给定的介质层数Medium和分层个数N1，对声速剖面重新插值，返回cP、cS和rho_k
 *
 * @param cp 输出的纵波速度数组
 * @param cs 输出的横波速度数组
 * @param rho_k 输出的密度数组
 * @param SSP 声速剖面结构体
 * @param Medium 介质层数
 * @param N1 分层个数
 */
void n2Linear(SSPStructure &SSP)
{
    int index = 0;
    for (size_t i = 0; i < SSP.NMedia; i++)
    {
        // 获取介质在SSP中的位置
        int ILoc = SSP.Loc(i);  // Fortran索引转C++索引
        int NPts = SSP.NPts(i); // 介质中的点数
        int N1 = SSP.N(i);      // 分层个数

        // 计算步长
        double h = (SSP.z(ILoc + NPts - 1) - SSP.z(ILoc)) / (N1 - 1);
        int Lay = 0; // 层索引

        // 遍历每个分层点
        for (int iz = 0; iz < N1; iz++)
        {
            int INDEX = index + iz;
            // 计算当前深度
            double z = SSP.z(ILoc) + iz * h;
            // 确保最后一个点的深度准确
            if (iz == N1 - 1)
            {
                z = SSP.z(ILoc + NPts - 1);
            }

            // 找到当前深度所在的层
            while (Lay < NPts - 1 && z > SSP.z(ILoc + Lay + 1))
            {
                Lay++;
            }

            // 计算插值参数
            int iSSP = ILoc + Lay;
            double dz = SSP.z(iSSP + 1) - SSP.z(iSSP);
            // 避免除零错误
            if (dz == 0.0)
            {
                dz = 1.0;
            }
            double R = (z - SSP.z(iSSP)) / dz;

            // 限制R在[0,1]范围内
            if (R < 0.0)
                R = 0.0;
            if (R > 1.0)
                R = 1.0;

            // P波速度计算 (N2线性插值)
            double alphaTop = SSP.alphaR(iSSP);
            double alphaBot = SSP.alphaR(iSSP + 1);
            // 避免除零错误
            if (alphaTop == 0.0)
                alphaTop = 1.0;
            if (alphaBot == 0.0)
                alphaBot = 1.0;

            std::complex<double> N2Top = 1.0 / (alphaTop * alphaTop);
            std::complex<double> N2Bot = 1.0 / (alphaBot * alphaBot);
            std::complex<double> N2Interp = (1.0 - R) * N2Top + R * N2Bot;
            // 避免开方错误
            if (real(N2Interp) <= 0.0)
            {
                SSP.cp_int(INDEX) = std::complex<double>(1500.0, 0.0); // 默认值
            }
            else
            {
                SSP.cp_int(INDEX) = 1.0 / sqrt(N2Interp);
            }

            // S波速度计算
            double betaTop = SSP.betaR(iSSP);
            double betaBot = SSP.betaR(iSSP + 1);
            if (betaTop != 0.0 && betaBot != 0.0)
            {
                N2Top = 1.0 / (betaTop * betaTop);
                N2Bot = 1.0 / (betaBot * betaBot);
                N2Interp = (1.0 - R) * N2Top + R * N2Bot;
                if (real(N2Interp) <= 0.0)
                {
                    SSP.cs_int(INDEX) = std::complex<double>(0.0, 0.0);
                }
                else
                {
                    SSP.cs_int(INDEX) = 1.0 / sqrt(N2Interp);
                }
            }
            else
            {
                SSP.cs_int(INDEX) = std::complex<double>(0.0, 0.0);
            }

            // 密度线性插值
            SSP.rho_int(INDEX) = (1.0 - R) * SSP.rho(iSSP) + R * SSP.rho(iSSP + 1);
        }
        index += N1;
    }
}

/**
 * @brief 根据给定的介质层数Medium和分层个数N1，对声速剖面重新插值，返回cP、cS和rho_k
 * Uses c-linear segments for P and S-wave speeds, rho-linear segments for density
 *
 * @param cp 输出的纵波速度数组
 * @param cs 输出的横波速度数组
 * @param rho_k 输出的密度数组
 * @param SSP 声速剖面结构体
 * @param Medium 介质层数
 * @param N1 分层个数
 */
void cLinear(SSPStructure &SSP)
{
    int index = 0;

    for (size_t i = 0; i < SSP.NMedia; i++)
    {
        // 获取介质在SSP中的位置
        int ILoc = SSP.Loc(i);  // Fortran索引转C++索引
        int NPts = SSP.NPts(i); // 介质中的点数
        int N1 = SSP.N(i);      // 分层个数

        // 计算步长
        double h = (SSP.z(ILoc + NPts - 1) - SSP.z(ILoc)) / (N1 - 1);
        int Lay = 0; // 层索引

        // 遍历每个分层点
        for (int iz = 0; iz < N1; iz++)
        {
            int INDEX = index + iz;
            // 计算当前深度
            double z = SSP.z(ILoc) + iz * h;
            // 确保最后一个点的深度准确
            if (iz == N1 - 1)
            {
                z = SSP.z(ILoc + NPts - 1);
            }

            // 找到当前深度所在的层
            while (Lay < NPts - 1 && z > SSP.z(ILoc + Lay + 1))
            {
                Lay++;
            }

            // 确保不会超出范围
            if (Lay >= NPts - 1)
            {
                Lay = NPts - 2;
            }

            // 计算插值参数
            int iSSP = ILoc + Lay;
            double dz = SSP.z(iSSP + 1) - SSP.z(iSSP);
            // 避免除零错误
            if (dz == 0.0)
            {
                dz = 1.0;
            }
            double R = (z - SSP.z(iSSP)) / dz;

            // 限制R在[0,1]范围内
            if (R < 0.0)
                R = 0.0;
            if (R > 1.0)
                R = 1.0;

            // P波速度计算 (c线性插值)
            SSP.cp_int(INDEX) = (1.0 - R) * SSP.cp(iSSP) + R * SSP.cp(iSSP + 1);

            // S波速度计算 (c线性插值)
            SSP.cs_int(INDEX) = (1.0 - R) * SSP.cs(iSSP) + R * SSP.cs(iSSP + 1);

            // 密度线性插值
            SSP.rho_int(INDEX) = (1.0 - R) * SSP.rho(iSSP) + R * SSP.rho(iSSP + 1);
        }
        index += N1;
    }
}

/**
 * @brief 根据给定的介质层数Medium和分层个数N1，对声速剖面重新插值，返回cP、cS和rho_k
 * Uses PCHIP segments for P, S-wave speeds and density rho
 *
 * @param cp 输出的纵波速度数组
 * @param cs 输出的横波速度数组
 * @param rho_k 输出的密度数组
 * @param SSP 声速剖面结构体
 * @param Medium 介质层数
 * @param N1 分层个数
 */
void cPCHIP(SSPStructure &SSP)
{
    int index = 0;

    for (size_t i = 0; i < SSP.NMedia; i++)
    {
        // 获取介质在SSP中的位置
        int ILoc = SSP.Loc(i);  // Fortran索引转C++索引
        int NPts = SSP.NPts(i); // 介质中的点数
        int N1 = SSP.N(i);      // 分层个数

        // 计算步长
        double h = (SSP.z(ILoc + NPts - 1) - SSP.z(ILoc)) / (N1 - 1);
        int Lay = 0; // 层索引

        // 遍历每个分层点
        for (int iz = 0; iz < N1; iz++)
        {
            int INDEX = index + iz;
            // 计算当前深度
            double z = SSP.z(ILoc) + iz * h;
            // 确保最后一个点的深度准确
            if (iz == N1 - 1)
            {
                z = SSP.z(ILoc + NPts - 1);
            }

            // 找到当前深度所在的层
            while (Lay < NPts - 1 && z > SSP.z(ILoc + Lay + 1))
            {
                Lay++;
            }

            // 确保不会超出范围
            if (Lay >= NPts - 1)
            {
                Lay = NPts - 2;
            }

            // 计算插值参数
            int iSSP = ILoc + Lay;
            double xt = z - SSP.z(iSSP);

            // P波速度计算 (PCHIP插值)
            SSP.cp_int(INDEX) = SSP.cpCoef(0, iSSP) +
                                (SSP.cpCoef(1, iSSP) +
                                 (SSP.cpCoef(2, iSSP) +
                                  SSP.cpCoef(3, iSSP) * xt) *
                                     xt) *
                                    xt;

            // S波速度计算 (PCHIP插值)
            SSP.cs_int(INDEX) = SSP.csCoef(0, iSSP) +
                                (SSP.csCoef(1, iSSP) +
                                 (SSP.csCoef(2, iSSP) +
                                  SSP.csCoef(3, iSSP) * xt) *
                                     xt) *
                                    xt;

            // 密度PCHIP插值
            SSP.rho_int(INDEX) = std::real(SSP.rhoCoef(0, iSSP) +
                                           (SSP.rhoCoef(1, iSSP) +
                                            (SSP.rhoCoef(2, iSSP) +
                                             SSP.rhoCoef(3, iSSP) * xt) *
                                                xt) *
                                               xt);
        }
        index += N1;
    }
}

/**
 * @brief 根据给定的介质层数Medium和分层个数N1，对声速剖面重新插值，返回cP、cS和rho_k
 * Uses cubic spline for P, S-wave speeds and density rho
 *
 * @param cp 输出的纵波速度数组
 * @param cs 输出的横波速度数组
 * @param rho_k 输出的密度数组
 * @param SSP 声速剖面结构体
 * @param Medium 介质层数
 * @param N1 分层个数
 */
void cCubic(SSPStructure &SSP)
{
    int index = 0;
    for (size_t i = 0; i < SSP.NMedia; i++)
    {
        // 获取介质在SSP中的位置
        int ILoc = SSP.Loc(i);  // Fortran索引转C++索引
        int NPts = SSP.NPts(i); // 介质中的点数
        int N1 = SSP.N(i);      // 分层个数

        // 计算步长
        double h = (SSP.z(ILoc + NPts - 1) - SSP.z(ILoc)) / (N1 - 1);
        int Lay = 0; // 层索引

        // 遍历每个分层点
        for (int iz = 0; iz < N1; iz++)
        {
            int INDEX = index + iz;
            // 计算当前深度
            double z = SSP.z(ILoc) + iz * h;
            // 确保最后一个点的深度准确
            if (iz == N1 - 1)
            {
                z = SSP.z(ILoc + NPts - 1);
            }

            // 找到当前深度所在的层
            while (Lay < NPts - 1 && z > SSP.z(ILoc + Lay + 1))
            {
                Lay++;
            }

            // 确保不会超出范围
            if (Lay >= NPts - 1)
            {
                Lay = NPts - 2;
            }

            // 计算插值参数
            int iSSP = ILoc + Lay;
            double xt = z - SSP.z(iSSP);

            // P波速度计算 (三次样条插值)
            SSP.cp_int(INDEX) = SSP.cpSpline(0, iSSP) +
                     (SSP.cpSpline(1, iSSP) +
                      (SSP.cpSpline(2, iSSP) +
                       SSP.cpSpline(3, iSSP) * xt) *
                          xt) *
                         xt;

            // S波速度计算 (三次样条插值)
            SSP.cs_int(INDEX) = SSP.csSpline(0, iSSP) +
                     (SSP.csSpline(1, iSSP) +
                      (SSP.csSpline(2, iSSP) +
                       SSP.csSpline(3, iSSP) * xt) *
                          xt) *
                         xt;

            // 密度三次样条插值
            SSP.rho_int(INDEX) = std::real(SSP.rhoSpline(0, iSSP) +
                                  (SSP.rhoSpline(1, iSSP) +
                                   (SSP.rhoSpline(2, iSSP) +
                                    SSP.rhoSpline(3, iSSP) * xt) *
                                       xt) *
                                      xt);
        }
        index += N1;
    }
}

/**
 * Munk profile
 *
 * Returns cs, cp, rho at depths i*h i = 1, N
 * Depths of interfaces
 *
 * @param cp Output complex vector for compressional wave speeds
 * @param cs Output complex vector for shear wave speeds
 * @param rho Output vector for densities
 * @param Medium Medium type (1 for ocean, 2 for fluid half-space, 9 for elastic layer)
 * @param N1 Number of points
 */
void Analytic(VectorXcd &cp, VectorXcd &cs, VectorXd &rho, int Medium, int N1)
{
    int N = N1 - 1;
    int i;
    double h, x, z;

    // Resize vectors to appropriate size
    cp.resize(N1);
    cs.resize(N1);
    rho.resize(N1);

    switch (Medium)
    {
    case 1: // THE OCEAN
        h = 5000.0 / N;
        for (i = 0; i < N1; i++)
        { // C++ uses 0-based indexing
            z = i * h;
            x = 2.0 * (z - 1300.0) / 1300.0;
            cp(i) = 1500.0 * (1.0 + eps * (x - 1.0 + exp(-x)));
            cs(i) = 0.0;
            rho(i) = 1.0;
        }
        break;

    case 2: // THE FLUID HALF-SPACE
        cp(0) = 1551.91;
        cs(0) = 0.0;
        rho(0) = 1.0e20;
        break;

    case 9: // AN ELASTIC LAYER
        h = 1000.0 / N;
        z = 5000.0;

        for (i = 0; i < N + 1; i++)
        { // C++ uses 0-based indexing
            cp(i) = 4700.0 + (z - 5000.0) / 10.0;
            cs(i) = 2000.0 + (z - 5000.0) / 10.0;
            cp(i) = 4700.0;
            cs(i) = 2000.0;
            rho(i) = 2.0;
            z = z + h;
        }
        break;
    }
}

void UpdateSSPLoss(double freq, double freq0, int Medium, SSPStructure &SSP)
{
    for (size_t i = 0; i < SSP.NMedia; i++)
    {
        int ILoc = SSP.Loc(Medium); // Fortran索引转C++索引
        int LocLen = SSP.NPts(Medium);
        for (size_t issp = 0; issp < SSP.NPts(Medium); issp++)
        {
            int iz = SSP.Loc(Medium) + issp;
            SSP.cp(iz) = CRCI(SSP.z(iz), SSP.alphaR(iz), SSP.alphaI(iz), freq, freq0,
                              SSP.AttenUnit, SSP.beta(Medium), SSP.ft(Medium));
            SSP.cp(iz) = CRCI(SSP.z(iz), SSP.alphaR(iz), SSP.alphaI(iz), freq, freq0,
                              SSP.AttenUnit, SSP.beta(Medium), SSP.ft(Medium));
            SSP.rho_int(iz) = SSP.rho(iz);
            SSP.cpSpline(1, iz) = SSP.cp(iz);
            SSP.csSpline(1, iz) = SSP.cs(iz);
            SSP.rhoSpline(1, iz) = SSP.rho(iz);
        }
        if (SSP.Type == SSP_Mode::MODE_P_cPCHIP)
        {
            PCHIP(SSP.z, SSP.cp, SSP.cpCoef, SSP.csWork, ILoc, LocLen);
            PCHIP(SSP.z, SSP.cs, SSP.csCoef, SSP.csWork, ILoc, LocLen);
            PCHIP(SSP.z, SSP.rho_int, SSP.rhoCoef, SSP.csWorkd, ILoc, LocLen);
        }
        if (SSP.Type == SSP_Mode::MODE_S_cCubic)
        {
            int IBCBeg = 0, IBCEnd = 0;
            CSpline(SSP.z, SSP.cpSpline, LocLen, IBCBeg, IBCEnd, LocLen, ILoc);
            CSpline(SSP.z, SSP.csSpline, LocLen, IBCBeg, IBCEnd, LocLen, ILoc);
            CSpline(SSP.z, SSP.rhoSpline, LocLen, IBCBeg, IBCEnd, LocLen, ILoc);
        }
    }
}

void UpdateHSLoss(double freq, double freq0, int Medium, SSPStructure &SSP, HSInfo &HS)
{
    double huge = 1e8;
    if (HS.BC == BC_Mode::MODE_A_Half_space)
    {
        HS.cp = CRCI(huge, HS.alphaR, HS.alphaI, freq, freq0,
                     SSP.AttenUnit, HS.beta, HS.ft);
        HS.cp = CRCI(huge, HS.betaR, HS.betaI, freq, freq0,
                     SSP.AttenUnit, HS.beta, HS.ft);
    }
}

void PCHIP(VectorXd &x, VectorXcd &y, MatrixXcd &PolyCoef, MatrixXcd &csWork, int istart, int nlen)
{
}
void PCHIP(VectorXd &x, VectorXd &y, MatrixXd &PolyCoef, MatrixXd &csWorkd, int istart, int nlen)
{
}
void CSpline(VectorXd &TAU, MatrixXcd &C, int N, int IBCBEG, int IBCEND, int NDIM, int istart)
{
}
void CSpline(VectorXd &TAU, MatrixXd &C, int N, int IBCBEG, int IBCEND, int NDIM, int istart)
{
}
