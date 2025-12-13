#include "BCImpedanceMod.h"

// 计算边界条件阻抗
void BCImpedance(size_t iprof, const double& x,  bool& isTop, complex<double> &f, complex<double> &g,
                 int &iPower, const bool& isComplex, KrakenMatrix &kramtrx, parameters& params, 
                 int& modeCount)
{
    int iTop = 0, iBot = 0;
    VectorXd yV = VectorXd::Zero(5);
    double mu;
    double rhoInside = 1.0;
    std::complex<double> gammaS, gammaP, gammaS2, gammaP2;
    std::complex<double> kx, kz, RCmplx, cInside(1500.0, 0.0);
    ReflectionCoef RInt;
    double omega = 2 * pi * params.freqinfo->freq;
    double omega2 = SQ(omega);
    double hFirstAcoustic = params.mesh.h(0);
    double hFirstAcoustic2 = SQ(hFirstAcoustic);
    HSInfo HS = params.HSTop[iprof];

    iPower = 0;

    // 获取边界内部的密度和声速
    if (isTop)
    {

        if (params.FirstAcoustic > -1)
        {
            iTop      = params.mesh.Loc( params.FirstAcoustic ) + params.mesh.N( params.FirstAcoustic );
            rhoInside = kramtrx.rho(iTop);
            cInside = std::sqrt(omega2 * hFirstAcoustic2) /
                      (2.0 + kramtrx.B1(0));
        }
    }
    else
    {
        if (params.LastAcoustic > -1)
        {
            iBot      = params.mesh.Loc( params.LastAcoustic ) + params.mesh.N( params.LastAcoustic );
            rhoInside = kramtrx.rho(iBot);
            cInside = std::sqrt(omega2 * hFirstAcoustic2 /
                      (2.0 + kramtrx.B1(kramtrx.B1.size() - 1)));
        }
        HS = params.HSBot[iprof];
    }

    // 根据边界条件类型返回阻抗
    if (HS.BC == BC_Mode::MODE_V_Vacuum)
    { // 真空边界
        f = 1.0;
        g = 0.0;
        yV(0) = std::real(f);
        yV(1) = std::real(g);
        yV(2) = 0.0;
        yV(3) = 0.0;
        yV(4) = 0.0;
    }
    else if (HS.BC == BC_Mode::MODE_R_Rigid)
    { // 刚性边界
        f = 0.0;
        g = 1.0;
        yV(0) = std::real(f);
        yV(1) = std::real(g);
        yV(2) = 0.0;
        yV(3) = 0.0;
        yV(4) = 0.0;
    }
    else if (HS.BC == BC_Mode::MODE_A_Half_space)
    { // 声弹性半空间
        if (real(HS.cs) > 0.0)
        {
            gammaS2 = x - omega2 / (real(HS.cs) * real(HS.cs));
            gammaP2 = x - omega2 / (real(HS.cp) * real(HS.cp));
            gammaS = std::sqrt(gammaS2);
            gammaP = std::sqrt(gammaP2);
            mu = HS.rho * real(HS.cs) * real(HS.cs);

            yV(0) = std::real((gammaS * gammaP - x) / mu);
            yV(1) = std::real(((gammaS2 + x) * (gammaS2 + x) - 4.0 * gammaS * gammaP * x) * mu);
            yV(2) = std::real(2.0 * gammaS * gammaP - gammaS2 - x);
            yV(3) = std::real(gammaP * (x - gammaS2));
            yV(4) = std::real(gammaS * (gammaS2 - x));

            f = omega2 * yV(3);
            g = yV(1);
            if (std::real(g) > 0.0)
            {
                modeCount += 1;
            }
        }
        else
        {
            gammaP = std::sqrt(std::complex<double>(x - omega2 / SQ(HS.cp)));
            f = gammaP;
            g = HS.rho;
            if (!isComplex)
            {
                f = std::real(f);
                g = std::real(g);
            }
        }
    }
    else if (HS.BC == BC_Mode::MODE_F_File)
    { // 表格化反射系数
        // 计算掠射角theta
        kx = std::sqrt(x);
        kz = std::sqrt(omega2 / (real(cInside) * real(cInside)) - x);
        RInt.theta = (180.0 / pi) * std::atan2(std::real(kz), std::real(kx)); // RadDeg的倒数

        // 计算R(ThetaInt)
        if (isTop)
        {
            InterpolateReflectionCoefficient(RInt, params.ReflectionCoef.RTop);
        }
        else
        {
            InterpolateReflectionCoefficient(RInt, params.ReflectionCoef.RBot);
        }

        // 将R(theta)转换为Robin边界条件中的(f,g)
        RCmplx = RInt.R * std::exp(I1D * RInt.phi);
        f = I1D * kz * (1.0 - RCmplx) / (rhoInside * (1.0 + RCmplx));
        g = 1.0;

        if (!isComplex)
        {
            f = 0.0;
            g = 1.0;
        }
    }
    else if (HS.BC == BC_Mode::MODE_P_Precomputed)
    { // 预计算反射系数
        // TODO InterpolateIRC

        if (!isComplex)
        {
            f = 0.0;
            g = 1.0;
            iPower = 0;
        }
    }

    // 顶部边界条件的符号与底部边界条件相反
    if (isTop)
    {
        g = -g;
    }

    // 穿过弹性层传播
    if (isTop)
    {
        if (params.FirstAcoustic > 1)
        { // 从顶部向下传播
            for (int im = 0; im < params.FirstAcoustic; ++im)
            {
                ElasticDN(x, yV, iPower, im, kramtrx, params);
            }

            f = omega2 * yV(3);
            g = yV(1);
        }
    }
    else
    {
        if (params.LastAcoustic < params.SSP[iprof].NMedia-1)
        { // 从底部向上传播
            for (int im = params.SSP[iprof].NMedia - 1; im > params.LastAcoustic; --im)
            {
                ElasticUP(x, yV, iPower, im, kramtrx, params);
            }

            f = omega2 * yV(3);
            g = yV(1);
        }
    }
}

// 向上传播通过弹性层
void ElasticUP(const double& x, VectorXd &yV, int &iPower, const int& Medium, KrakenMatrix &kramtrx, parameters& params)
{
    VectorXd xV(5), zV(5);

    double h = params.mesh.h(Medium);
    // 第一步使用欧拉法
    double two_x = 2.0 * x;
    double two_h = 2.0 * h;
    double four_h_x = 4.0 * h * x;
    int j = params.mesh.Loc(Medium) + params.mesh.N(Medium); // C++是0-based
    double xB3 = x * kramtrx.B3(j) - kramtrx.B1(j);

    zV(0) = yV(0) - 0.5 * (kramtrx.B1(j) * yV(3) - kramtrx.B2(j) * yV(4));
    zV(1) = yV(1) - 0.5 * (-kramtrx.rho(j) * yV(3) - xB3 * yV(4));
    zV(2) = yV(2) - 0.5 * (two_h * yV(3) + kramtrx.B4(j) * yV(4));
    zV(3) = yV(3) - 0.5 * (xB3 * yV(0) + kramtrx.B2(j) * yV(1) - two_x * kramtrx.B4(j) * yV(2));
    zV(4) = yV(4) - 0.5 * (kramtrx.rho(j) * yV(0) - kramtrx.B1(j) * yV(1) - four_h_x * yV(2));

    // 改进的中点法
    for (size_t ii = params.mesh.N(Medium) - 1; ii >= 0; --ii)
    {
        j--;

        // 保存前一步的值
        for (int k = 0; k < 5; ++k)
        {
            xV(k) = yV(k);
        }

        // 更新当前值
        for (int k = 0; k < 5; ++k)
        {
            yV(k) = zV(k);
        }

        xB3 = x * kramtrx.B3(j) - kramtrx.B1(j);

        zV(0) = xV(0) - (kramtrx.B1(j) * yV(3) - kramtrx.B2(j) * yV(4));
        zV(1) = xV(1) - (-kramtrx.rho(j) * yV(3) - xB3 * yV(4));
        zV(2) = xV(2) - (two_h * yV(3) + kramtrx.B4(j) * yV(4));
        zV(3) = xV(3) - (xB3 * yV(0) + kramtrx.B2(j) * yV(1) - two_x * kramtrx.B4(j) * yV(2));
        zV(4) = xV(4) - (kramtrx.rho(j) * yV(0) - kramtrx.B1(j) * yV(1) - four_h_x * yV(2));

        // 必要时进行缩放
        if (ii != 0)
        {
            if (std::abs(zV(1)) < BCIFloor)
            {
                for (int k = 0; k < 5; ++k)
                {
                    zV(k) *= BCIRoof;
                    yV(k) *= BCIRoof;
                }
                iPower -= BCIiPowerR;
            }

            if (std::abs(zV(1)) > BCIRoof)
            {
                for (int k = 0; k < 5; ++k)
                {
                    zV(k) *= BCIFloor;
                    yV(k) *= BCIFloor;
                }
                iPower -= BCIiPowerF;
            }
        }
    }

    // 在终点应用标准滤波器
    for (int k = 0; k < 5; ++k)
    {
        yV(k) = (xV(k) + 2.0 * yV(k) + zV(k)) / 4.0;
    }
}

// 向下传播通过弹性层
void ElasticDN(const double x, VectorXd &yV, int &iPower, const int& Medium, KrakenMatrix &kramtrx, parameters& params)
{
    VectorXd xV(5), zV(5);

    // 第一步使用欧拉法
    double two_x = 2.0 * x;
    double two_h = 2.0 * params.mesh.h(Medium);
    double four_h_x = 4.0 * params.mesh.h(Medium) * x;
    int j = params.mesh.Loc(Medium); // C++是0-based
    double xB3 = x * kramtrx.B3(j) - kramtrx.B1(j);

    zV(0) = yV(0) + 0.5 * (kramtrx.B1(j) * yV(3) - kramtrx.B2(j) * yV(4));
    zV(1) = yV(1) + 0.5 * (-kramtrx.rho(j) * yV(3) - xB3 * yV(4));
    zV(2) = yV(2) + 0.5 * (two_h * yV(3) + kramtrx.B4(j) * yV(4));
    zV(3) = yV(3) + 0.5 * (xB3 * yV(0) + kramtrx.B2(j) * yV(1) - two_x * kramtrx.B4(j) * yV(2));
    zV(4) = yV(4) + 0.5 * (kramtrx.rho(j) * yV(0) - kramtrx.B1(j) * yV(1) - four_h_x * yV(2));

    // 改进的中点法
    for (size_t ii = 0; ii < params.mesh.N(Medium); ++ii)
    {
        j++;

        // 保存前一步的值
        for (int k = 0; k < 5; ++k)
        {
            xV(k) = yV(k);
        }

        // 更新当前值
        for (int k = 0; k < 5; ++k)
        {
            yV(k) = zV(k);
        }

        xB3 = x * kramtrx.B3(j) - kramtrx.B1(j);

        zV(0) = xV(0) + (kramtrx.B1(j) * yV(3) - kramtrx.B2(j) * yV(4));
        zV(1) = xV(1) + (-kramtrx.rho(j) * yV(3) - xB3 * yV(4));
        zV(2) = xV(2) + (two_h * yV(3) + kramtrx.B4(j) * yV(4));
        zV(3) = xV(3) + (xB3 * yV(0) + kramtrx.B2(j) * yV(1) - two_x * kramtrx.B4(j) * yV(2));
        zV(4) = xV(4) + (kramtrx.rho(j) * yV(0) - kramtrx.B1(j) * yV(1) - four_h_x * yV(2));

        // 必要时进行缩放
        if (ii != params.mesh.N(Medium) - 1)
        {
            if (std::abs(zV(1)) < BCIFloor)
            {
                for (int k = 0; k < 5; ++k)
                {
                    zV(k) *= BCIRoof;
                    yV(k) *= BCIRoof;
                }
                iPower -= BCIiPowerR;
            }

            if (std::abs(zV(1)) > BCIRoof)
            {
                for (int k = 0; k < 5; ++k)
                {
                    zV(k) *= BCIFloor;
                    yV(k) *= BCIFloor;
                }
                iPower -= BCIiPowerF;
            }
        }
    }

    // 在终点应用标准滤波器
    for (int k = 0; k < 5; ++k)
    {
        yV(k) = (xV(k) + 2.0 * yV(k) + zV(k)) / 4.0;
    }
}