#include "sspMod.h"

void EvaluateSSP(SSPStructure &SSP, SSP_Mode &ssptype)
{
    // double cxx, cyy, cxy, cxz, cyz;

    switch (ssptype)
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
        std::cerr << "Unknown SSP type in EvaluateSSP, use default type *cLinear*" << std::endl;
        cLinear(SSP);
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
    // 计算步长
    double h = (SSP.z(SSP.z.size() - 1) - SSP.z(0)) / SSP.N;
    int Lay = 0; // 层索引
    SSP.cp_int.resize(SSP.N + 1);
    SSP.cs_int.resize(SSP.N + 1);
    SSP.rho_int.resize(SSP.N + 1);

    // 遍历每个分层点
    for (int iz = 0; iz < SSP.N + 1; iz++)
    {
        // 计算当前深度
        double z = SSP.z(0) + iz * h;
        // 确保最后一个点的深度准确
        if (iz == SSP.N + 1 - 1)
        {
            z = SSP.z(SSP.z.size() - 1);
        }

        // 找到当前深度所在的层
        while (Lay < SSP.NPts - 1 && z > SSP.z(Lay + 1))
        {
            Lay++;
        }

        // 计算插值参数
        double dz = SSP.z(Lay + 1) - SSP.z(Lay);
        // 避免除零错误
        if (dz == 0.0)
        {
            dz = 1.0;
        }
        double R = (z - SSP.z(Lay)) / dz;

        // 限制R在[0,1]范围内
        if (R < 0.0)
            R = 0.0;
        if (R > 1.0)
            R = 1.0;

        // P波速度计算 (N2线性插值)
        complex<double> alphaTop = SSP.cp(Lay);
        complex<double> alphaBot = SSP.cp(Lay + 1);
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
            SSP.cp_int(iz) = std::complex<double>(1500.0, 0.0); // 默认值
        }
        else
        {
            SSP.cp_int(iz) = 1.0 / sqrt(N2Interp);
        }

        // S波速度计算
        complex<double> betaTop = SSP.cs(Lay);
        complex<double> betaBot = SSP.cs(Lay + 1);
        if (betaTop != 0.0 && betaBot != 0.0)
        {
            N2Top = 1.0 / (betaTop * betaTop);
            N2Bot = 1.0 / (betaBot * betaBot);
            N2Interp = (1.0 - R) * N2Top + R * N2Bot;
            if (real(N2Interp) <= 0.0)
            {
                SSP.cs_int(iz) = std::complex<double>(0.0, 0.0);
            }
            else
            {
                SSP.cs_int(iz) = 1.0 / sqrt(N2Interp);
            }
        }
        else
        {
            SSP.cs_int(iz) = std::complex<double>(0.0, 0.0);
        }

        // 密度线性插值
        SSP.rho_int(iz) = (1.0 - R) * SSP.rho(Lay) + R * SSP.rho(Lay + 1);
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
    // 计算步长
    double h = (SSP.z(SSP.z.size() - 1) - SSP.z(0)) / SSP.N;
    int Lay = 0; // 层索引
    SSP.cp_int.resize(SSP.N + 1);
    SSP.cs_int.resize(SSP.N + 1 + 1);
    SSP.rho_int.resize(SSP.N + 1);

    // 遍历每个分层点
    for (int iz = 0; iz < SSP.N + 1; iz++)
    {
        // 计算当前深度
        double z = SSP.z(0) + iz * h;
        // 确保最后一个点的深度准确
        if (iz == SSP.N + 1 - 1)
        {
            z = SSP.z(SSP.z.size() - 1);
        }

        // 找到当前深度所在的层
        while (Lay < SSP.NPts - 1 && z > SSP.z(Lay + 1))
        {
            Lay++;
        }

        // 确保不会超出范围
        if (Lay >= SSP.NPts - 1)
        {
            Lay = SSP.NPts - 2;
        }

        // 计算插值参数
        double dz = SSP.z(Lay + 1) - SSP.z(Lay);
        // 避免除零错误
        if (dz == 0.0)
        {
            dz = 1.0;
        }
        double R = (z - SSP.z(Lay)) / dz;

        // 限制R在[0,1]范围内
        if (R < 0.0)
            R = 0.0;
        if (R > 1.0)
            R = 1.0;

        // P波速度计算 (c线性插值)
        SSP.cp_int(iz) = (1.0 - R) * SSP.cp(Lay) + R * SSP.cp(Lay + 1);

        // S波速度计算 (c线性插值)
        SSP.cs_int(iz) = (1.0 - R) * SSP.cs(Lay) + R * SSP.cs(Lay + 1);

        // 密度线性插值
        SSP.rho_int(iz) = (1.0 - R) * SSP.rho(Lay) + R * SSP.rho(Lay + 1);
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
    // 计算步长
    double h = (SSP.z(SSP.z.size() - 1) - SSP.z(0)) / SSP.N;
    int Lay = 0; // 层索引
    SSP.cp_int.resize(SSP.N + 1);
    SSP.cs_int.resize(SSP.N + 1);
    SSP.rho_int.resize(SSP.N + 1);

    // 遍历每个分层点
    for (int iz = 0; iz < SSP.N + 1; iz++)
    {
        // 计算当前深度
        double z = SSP.z(0) + iz * h;
        // 确保最后一个点的深度准确
        if (iz == SSP.N + 1 - 1)
        {
            z = SSP.z(SSP.z.size() - 1);
        }

        // 找到当前深度所在的层
        while (Lay < SSP.NPts - 1 && z > SSP.z(Lay + 1))
        {
            Lay++;
        }

        // 确保不会超出范围
        if (Lay >= SSP.NPts - 1)
        {
            Lay = SSP.NPts - 2;
        }

        // 计算插值参数
        double xt = z - SSP.z(Lay);

        // P波速度计算 (PCHIP插值)
        SSP.cp_int(iz) = SSP.cpCoef(0, Lay) +
                         (SSP.cpCoef(1, Lay) +
                          (SSP.cpCoef(2, Lay) +
                           SSP.cpCoef(3, Lay) * xt) *
                              xt) *
                             xt;

        // S波速度计算 (PCHIP插值)
        SSP.cs_int(iz) = SSP.csCoef(0, Lay) +
                         (SSP.csCoef(1, Lay) +
                          (SSP.csCoef(2, Lay) +
                           SSP.csCoef(3, Lay) * xt) *
                              xt) *
                             xt;

        // 密度PCHIP插值
        SSP.rho_int(iz) = std::real(SSP.rhoCoef(0, Lay) +
                                    (SSP.rhoCoef(1, Lay) +
                                     (SSP.rhoCoef(2, Lay) +
                                      SSP.rhoCoef(3, Lay) * xt) *
                                         xt) *
                                        xt);
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
    // 计算步长
    double h = (SSP.z(SSP.z.size() - 1) - SSP.z(0)) / SSP.N;
    int Lay = 0; // 层索引
    std::complex<double> cp_cmplx, cs_cmplx, rho_cmplx, cpz_cmplx, cpzz_cmplx;
    SSP.cp_int.resize(SSP.N + 1);
    SSP.cs_int.resize(SSP.N + 1);
    SSP.rho_int.resize(SSP.N + 1);

    // 遍历每个分层点
    for (int iz = 0; iz < SSP.N + 1; iz++)
    {
        // 计算当前深度
        double z = SSP.z(0) + iz * h;
        // 确保最后一个点的深度准确
        if (iz == SSP.N + 1 - 1)
        {
            z = SSP.z(SSP.z.size() - 1);
        }

        // 找到当前深度所在的层
        while (Lay < SSP.NPts - 1 && z > SSP.z(Lay + 1))
        {
            Lay++;
        }

        // 确保不会超出范围
        if (Lay >= SSP.NPts - 1)
        {
            Lay = SSP.NPts - 2;
        }

        // 计算插值参数
        double xt = z - SSP.z(Lay);

        SplineALL(SSP.cpSpline, Lay, xt, cp_cmplx, cpz_cmplx, cpzz_cmplx);
        SplineALL(SSP.csSpline, Lay, xt, cs_cmplx, cpz_cmplx, cpzz_cmplx);
        SplineALL(SSP.rhoSpline, Lay, xt, rho_cmplx, cpz_cmplx, cpzz_cmplx);
        SSP.cp_int(iz) = cp_cmplx;
        SSP.cs_int(iz) = cs_cmplx;
        SSP.rho_int(iz) = real(rho_cmplx);

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

void UpdateSSPLoss(double &freq, double &freq0, int &NMedia, SSP_Mode &SSPType, Atten_Mode &AttenUnit, SSPStructure *SSPList)
{
    for (size_t i = 0; i < NMedia; i++)
    {
        SSPStructure &SSP = SSPList[i];
        for (size_t iz = 0; iz < SSP.NPts; iz++)
        {
            SSP.cp(iz) = CRCI(SSP.z(iz), SSP.alphaR(iz), SSP.alphaI(iz), freq, freq0,
                              AttenUnit, SSP.beta, SSP.ft);
            SSP.cs(iz) = CRCI(SSP.z(iz), SSP.betaR(iz), SSP.betaI(iz), freq, freq0,
                              AttenUnit, SSP.beta, SSP.ft);
            if (SSPType == SSP_Mode::MODE_S_cCubic)
            {
                SSP.cpSpline(0, iz) = SSP.cp(iz);
                SSP.csSpline(0, iz) = SSP.cs(iz);
                SSP.rhoSpline(0, iz) = SSP.rho(iz);
            }
        }
        if (SSPType == SSP_Mode::MODE_P_cPCHIP)
        {
            PCHIP(SSP.z, SSP.cp, SSP.NPts, SSP.cpCoef, SSP.csWork);
            PCHIP(SSP.z, SSP.cs, SSP.NPts, SSP.csCoef, SSP.csWork);
            VectorXcd rhoc = SSP.rho.cast<std::complex<double>>();
            PCHIP(SSP.z, rhoc, SSP.NPts, SSP.rhoCoef, SSP.csWork);
        }
        if (SSPType == SSP_Mode::MODE_S_cCubic)
        {
            int IBCBeg = 0, IBCEnd = 0;
            CSpline(SSP.z, SSP.cpSpline, SSP.NPts, IBCBeg, IBCEnd, SSP.NPts);
            CSpline(SSP.z, SSP.csSpline, SSP.NPts, IBCBeg, IBCEnd, SSP.NPts);
            CSpline(SSP.z, SSP.rhoSpline, SSP.NPts, IBCBeg, IBCEnd, SSP.NPts);
        }
    }
}

void UpdateHSLoss(double &freq, double &freq0, int &Medium, Atten_Mode &AttenUnit, HSInfo &HSTop, HSInfo &HSBot)
{
    double huge = 1e8;
    if (HSTop.BC == BC_Mode::MODE_A_Half_space)
    {
        HSTop.cp = CRCI(huge, HSTop.alphaR, HSTop.alphaI, freq, freq0,
                        AttenUnit, HSTop.beta, HSTop.ft);
        HSTop.cs = CRCI(huge, HSTop.betaR, HSTop.betaI, freq, freq0,
                        AttenUnit, HSTop.beta, HSTop.ft);
    }
    if (HSBot.BC == BC_Mode::MODE_A_Half_space)
    {
        HSBot.cp = CRCI(huge, HSBot.alphaR, HSBot.alphaI, freq, freq0,
                        AttenUnit, HSBot.beta, HSBot.ft);
        HSBot.cs = CRCI(huge, HSBot.betaR, HSBot.betaI, freq, freq0,
                        AttenUnit, HSBot.beta, HSBot.ft);
    }
}
