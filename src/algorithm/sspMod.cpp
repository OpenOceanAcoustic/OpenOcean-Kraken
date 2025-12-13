#include "sspMod.h"

void EvaluateSSP(SSPStructure &SSP, SSP_Mode &ssptype,int iMedium)
{
    // double cxx, cyy, cxy, cxz, cyz;
    // 输入校验
    assert(iMedium >= 0 && iMedium < SSP.NMedia);
    switch (ssptype)
    {
    case SSP_Mode::MODE_N_n2Linear:
        n2Linear(SSP,iMedium);
        break;
    case SSP_Mode::MODE_C_cLinear:
        cLinear(SSP,iMedium);
        break;
    case SSP_Mode::MODE_P_cPCHIP:
        cPCHIP(SSP,iMedium);
        break;
    case SSP_Mode::MODE_S_cCubic:
        cCubic(SSP,iMedium);
        break;
    // case SSP_Mode::MODE_A_Analytic:
    //     Analytic(cp, cs, rho_k, Medium, N1); // 需要实现Analytic函数
    //     break;
    default:
        // 报错
        std::cerr << "Unknown SSP type in EvaluateSSP, use default type *cLinear*" << std::endl;
        cLinear(SSP,iMedium);
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
void n2Linear(SSPStructure &SSP,int iMedium)
{   
    // 确定介质位置
    int offset_start = SSP.get_media_start(iMedium);
    int offset_end = SSP.get_media_end(iMedium);
    int current_medium_Nmesh = SSP.get_media_Nmesh(iMedium); // 当前介质网格点数
    int current_medium_size = SSP.get_media_size(iMedium); // 当前介质点数
    // 计算当前介质在全局插值数组中的起始位置
    int global_offset = SSP.get_global_interp_offset(iMedium);
    // 计算步长
    double h = (SSP.z(offset_end) - SSP.z(offset_start)) / current_medium_Nmesh;
    int Lay = 0; // 层索引
    // SSP.cp_int.resize(SSP.N + 1);
    // SSP.cs_int.resize(SSP.N + 1);
    // SSP.rho_int.resize(SSP.N + 1);

    // 遍历每个分层点
    for (int local_iz  = 0; local_iz  <= current_medium_Nmesh + 1; local_iz ++)
    {
        // 计算当前深度
        double z = SSP.z(offset_start) + local_iz  * h;
        // 确保最后一个点的深度准确
        if (local_iz  == current_medium_Nmesh)
        {
            z = SSP.z(offset_end);
        }

        // 找到当前深度所在的层
        while (Lay < current_medium_size - 1 && z > SSP.z(offset_start+Lay + 1))
        {
            Lay++;
        }

        // 计算插值参数 R ∈ [0,1]
        double z_top = SSP.z[offset_start + Lay];
        double z_bot = SSP.z[offset_start + Lay + 1];
        double dz = z_bot - z_top;
        if (std::abs(dz) < 1e-12) dz = 1.0; // 避免除零
        double R = (z - z_top) / dz;
        R = std::max(0.0, std::min(1.0, R)); // clamp to [0,1]

        // P波速度计算 (N2线性插值)
        std::complex<double> alphaTop =SSP.cp[offset_start + Lay];
        std::complex<double> alphaBot = SSP.cp[offset_start + Lay + 1];
        // 防止除零（物理上 cp 不应为 0，但做安全处理）
        if (std::abs(alphaTop) < 1e-12) alphaTop = std::complex<double>(1500.0, 0.0);
        if (std::abs(alphaBot) < 1e-12) alphaBot = std::complex<double>(1500.0, 0.0);

        std::complex<double> N2Top = 1.0 / (alphaTop * alphaTop);
        std::complex<double> N2Bot = 1.0 / (alphaBot * alphaBot);
        std::complex<double> N2Interp = (1.0 - R) * N2Top + R * N2Bot;

         int global_idx = global_offset + local_iz;
        // 避免开方错误
       if (std::real(N2Interp) <= 0.0) {
            SSP.cp_int[global_idx] = std::complex<double>(1500.0, 0.0); // 默认声速
        } else {
            SSP.cp_int[global_idx] = 1.0 / std::sqrt(N2Interp);
        }


        // S波速度计算
         std::complex<double> betaTop = SSP.cs[offset_start + Lay];
        std::complex<double> betaBot = SSP.cs[offset_start + Lay + 1];

         if (std::abs(betaTop) > 1e-12 && std::abs(betaBot) > 1e-12) {
            N2Top = 1.0 / (betaTop * betaTop);
            N2Bot = 1.0 / (betaBot * betaBot);
            N2Interp = (1.0 - R) * N2Top + R * N2Bot;
            if (std::real(N2Interp) <= 0.0) {
                SSP.cs_int[global_idx] = std::complex<double>(0.0, 0.0);
            } else {
                SSP.cs_int[global_idx] = 1.0 / std::sqrt(N2Interp);
            }
        } else {
            SSP.cs_int[global_idx] = std::complex<double>(0.0, 0.0);
        }

        // --- 密度：线性插值 ---
        double rhoTop = SSP.rho[offset_start + Lay];
        double rhoBot = SSP.rho[offset_start + Lay + 1];
        SSP.rho_int[global_idx] = (1.0 - R) * rhoTop + R * rhoBot;
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
void cLinear(SSPStructure &SSP, int iMedium)
{

    // 获取当前介质在原始数据中的范围（注意：get_media_end 返回的是半开区间 end）
    int orig_start = SSP.get_media_start(iMedium);
    int orig_end   = SSP.get_media_end(iMedium);      // = orig_start + NPts[iMedium]
    int orig_size  = SSP.get_media_size(iMedium);     // = NPts[iMedium]
    int Nmesh      = SSP.get_media_Nmesh(iMedium);    // 插值段数

    // 计算当前介质在插值输出数组中的起始位置
    int interp_start = SSP.get_global_interp_offset(iMedium);

    // 确保插值数组已分配足够空间（可选：也可由调用者保证）
    // int total_interp_points = 0;
    // for (int i = 0; i < SSP.NMedia; ++i) {
    //     total_interp_points += SSP.NMesh[i] + 1;
    // }
    // if (static_cast<int>(SSP.cp_int.size()) != total_interp_points) {
    //     SSP.cp_int.resize(total_interp_points);
    //     SSP.cs_int.resize(total_interp_points);
    //     SSP.rho_int.resize(total_interp_points);
    // }

    // 当前介质的深度范围
    double z_top = SSP.z[orig_start];
    double z_bot = SSP.z[orig_end - 1]; // 因为 orig_end 是半开区间
    double h = (z_bot - z_top) / static_cast<double>(Nmesh);

    int Lay = 0; // 在当前介质内的局部层索引（0 ~ orig_size-2）

    // 遍历插值点：共 Nmesh + 1 个点（对应 iz = 0 到 Nmesh）
    for (int local_iz = 0; local_iz <= Nmesh; ++local_iz)
    {
        // 计算当前深度 z
        double z = z_top + local_iz * h;
        if (local_iz == Nmesh) {
            z = z_bot; // 精确保底，避免浮点误差导致越界
        }

        // 在当前介质内查找 z 所在的原始子层
        while (Lay < orig_size - 1 && z > SSP.z[orig_start + Lay + 1]) {
            Lay++;
        }

        // 防止 Lay 越界（理论上不会发生，但安全起见）
        if (Lay >= orig_size - 1) Lay = orig_size - 2;
        if (Lay < 0) Lay = 0;

        // 插值参数 R ∈ [0, 1]
        double zL = SSP.z[orig_start + Lay];
        double zR = SSP.z[orig_start + Lay + 1];
        double dz = zR - zL;
        if (std::abs(dz) < 1e-14) dz = 1.0; // 避免除零

        double R = (z - zL) / dz;
        R = std::max(0.0, std::min(1.0, R)); // clamp to [0,1]

        // 全局输出索引
        int global_idx = interp_start + local_iz;

        // --- P 波速度：c-linear ---
        SSP.cp_int[global_idx] = (1.0 - R) * SSP.cp[orig_start + Lay] + R * SSP.cp[orig_start + Lay + 1];

        // --- S 波速度：c-linear ---
        SSP.cs_int[global_idx] = (1.0 - R) * SSP.cs[orig_start + Lay] + R * SSP.cs[orig_start + Lay + 1];

        // --- 密度：线性 ---
        SSP.rho_int[global_idx] = (1.0 - R) * SSP.rho[orig_start + Lay] + R * SSP.rho[orig_start + Lay + 1];
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
void cPCHIP(SSPStructure &SSP, int iMedium)
{
    assert(iMedium >= 0 && iMedium < SSP.NMedia);

    // 获取当前介质在原始数据中的范围
    int orig_start = SSP.get_media_start(iMedium);
    int orig_end   = SSP.get_media_end(iMedium);      // 半开区间 [start, end)
    int orig_size  = SSP.get_media_size(iMedium);     // = NPts[iMedium]
    int Nmesh      = SSP.get_media_Nmesh(iMedium);    // 插值段数

    // 计算当前介质在插值输出数组中的起始位置
    int interp_start = SSP.get_global_interp_offset(iMedium);

    // 预分配插值数组（如果尚未分配）
    // int total_interp = 0;
    // for (int i = 0; i < SSP.NMedia; ++i) {
    //     total_interp += SSP.NMesh[i] + 1;
    // }
    // if (static_cast<int>(SSP.cp_int.size()) != total_interp) {
    //     SSP.cp_int.resize(total_interp);
    //     SSP.cs_int.resize(total_interp);
    //     SSP.rho_int.resize(total_interp);
    // }

    // 当前介质的深度范围
    double z_top = SSP.z[orig_start];
    double z_bot = SSP.z[orig_end - 1];
    double h = (z_bot - z_top) / static_cast<double>(Nmesh);

    int Lay = 0; // 在当前介质内的局部层索引（0 ~ orig_size-2）

    // 遍历插值点：共 Nmesh + 1 个
    for (int local_iz = 0; local_iz <= Nmesh; ++local_iz)
    {
        double z = z_top + local_iz * h;
        if (local_iz == Nmesh) {
            z = z_bot; // 精确保底
        }

        // 在当前介质内查找 z 所在的原始子层
        while (Lay < orig_size - 1 && z > SSP.z[orig_start + Lay + 1]) {
            Lay++;
        }
        if (Lay >= orig_size - 1) Lay = orig_size - 2;
        if (Lay < 0) Lay = 0;

        // 相对于当前子层起点的偏移
        double xt = z - SSP.z[orig_start + Lay];

        // 全局系数列索引 = orig_start + Lay
        int coef_col = orig_start + Lay;

        // 全局输出索引
        int global_idx = interp_start + local_iz;

        // --- P 波速度：PCHIP ---
        std::complex<double> cp_val =
            SSP.cpCoef(0, coef_col) +
            (SSP.cpCoef(1, coef_col) +
             (SSP.cpCoef(2, coef_col) +
              SSP.cpCoef(3, coef_col) * xt) * xt) * xt;
        SSP.cp_int[global_idx] = cp_val;

        // --- S 波速度：PCHIP ---
        std::complex<double> cs_val =
            SSP.csCoef(0, coef_col) +
            (SSP.csCoef(1, coef_col) +
             (SSP.csCoef(2, coef_col) +
              SSP.csCoef(3, coef_col) * xt) * xt) * xt;
        SSP.cs_int[global_idx] = cs_val;

        // --- 密度：PCHIP（取实部）---
        std::complex<double> rho_val =
            SSP.rhoCoef(0, coef_col) +
            (SSP.rhoCoef(1, coef_col) +
             (SSP.rhoCoef(2, coef_col) +
              SSP.rhoCoef(3, coef_col) * xt) * xt) * xt;
        SSP.rho_int[global_idx] = std::real(rho_val);
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
void cCubic(SSPStructure &SSP, int iMedium)
{
    

    // 获取当前介质在原始数据中的范围
    int orig_start = SSP.get_media_start(iMedium);
    int orig_end   = SSP.get_media_end(iMedium);      // 半开区间 [start, end)
    int orig_size  = SSP.get_media_size(iMedium);     // = NPts[iMedium]
    int Nmesh      = SSP.get_media_Nmesh(iMedium);    // 插值段数

    // 计算当前介质在插值输出数组中的起始位置
    int interp_start = SSP.get_global_interp_offset(iMedium);

    // // 预分配插值数组（如果尚未分配）
    // int total_interp = 0;
    // for (int i = 0; i < SSP.NMedia; ++i) {
    //     total_interp += SSP.NMesh[i] + 1;
    // }
    // if (static_cast<int>(SSP.cp_int.size()) != total_interp) {
    //     SSP.cp_int.resize(total_interp);
    //     SSP.cs_int.resize(total_interp);
    //     SSP.rho_int.resize(total_interp);
    // }

    // 当前介质的深度范围
    double z_top = SSP.z[orig_start];
    double z_bot = SSP.z[orig_end - 1];
    double h_step = (z_bot - z_top) / static_cast<double>(Nmesh);

    int Lay = 0; // 在当前介质内的局部层索引（0 ~ orig_size-2）

    // 遍历插值点：共 Nmesh + 1 个
    for (int local_iz = 0; local_iz <= Nmesh; ++local_iz)
    {
        double z = z_top + local_iz * h_step;
        if (local_iz == Nmesh) {
            z = z_bot; // 精确保底
        }

        // 在当前介质内查找 z 所在的原始子层
        while (Lay < orig_size - 1 && z > SSP.z[orig_start + Lay + 1]) {
            Lay++;
        }
        if (Lay >= orig_size - 1) Lay = orig_size - 2;
        if (Lay < 0) Lay = 0;

        // 相对于当前子层起点的偏移（即 H）
        double H = z - SSP.z[orig_start + Lay];

        // 全局样条系数列索引
        int coef_col = orig_start + Lay;

        // 全局输出索引
        int global_idx = interp_start + local_iz;

        // --- P 波速度：三次样条 ---
        std::complex<double> cp_val, cpz, cpzz;
        SplineALL(SSP.cpSpline, coef_col, H, cp_val, cpz, cpzz);
        SSP.cp_int[global_idx] = cp_val;

        // --- S 波速度：三次样条 ---
        std::complex<double> cs_val, csz, cszz;
        SplineALL(SSP.csSpline, coef_col, H, cs_val, csz, cszz);
        SSP.cs_int[global_idx] = cs_val;

        // --- 密度：三次样条（取实部）---
        std::complex<double> rho_val, rhoz, rhozz;
        SplineALL(SSP.rhoSpline, coef_col, H, rho_val, rhoz, rhozz);
        SSP.rho_int[global_idx] = std::real(rho_val);
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

void UpdateSSPLoss(double freq, double freq0,
                   int NMedia, SSP_Mode SSPType, Atten_Mode AttenUnit,
                   SSPStructure& ssp)
{
    for (int iMedium = 0; iMedium < NMedia; ++iMedium)
    {
        const int start = ssp.get_media_start(iMedium);
        const int npoints = ssp.get_media_size(iMedium);

        // Step 1: 更新 cp, cs, rhoT 在该层的所有点
        for (int idx = 0; idx < npoints; ++idx)
        {
            const int iz = start + idx; // 全局索引

            ssp.cp(iz) = CRCI(ssp.z(iz), ssp.alphaR(iz), ssp.alphaI(iz),
                              freq, freq0, AttenUnit,
                              ssp.beta[iMedium], ssp.ft[iMedium]);

            ssp.cs(iz) = CRCI(ssp.z(iz), ssp.betaR(iz), ssp.betaI(iz),
                              freq, freq0, AttenUnit,
                              ssp.beta[iMedium], ssp.ft[iMedium]);


        }

        // Step 2: 如果使用 PCHIP 插值
        if (SSPType == SSP_Mode::MODE_P_cPCHIP)
        {
            // 提取当前层的子向量
            Eigen::VectorXd z_seg = ssp.z.segment(start, npoints);
            Eigen::VectorXcd cp_seg = ssp.cp.segment(start, npoints);
            Eigen::VectorXcd cs_seg = ssp.cs.segment(start, npoints);
            Eigen::VectorXcd rho_seg = ssp.rho.segment(start, npoints);
            
            // 创建局部系数矩阵（4 x npoints）
            Eigen::MatrixXcd cpCoef_local(4, npoints);
            Eigen::MatrixXcd csCoef_local(4, npoints);
            Eigen::MatrixXcd rhoCoef_local(4, npoints);
            Eigen::MatrixXcd work_local(2, npoints); // 假设 csWork 至少需要 2 行

            // 调用 PCHIP（传完整局部矩阵，非 block）
            PCHIP(z_seg, cp_seg, npoints, cpCoef_local, work_local);
            PCHIP(z_seg, cs_seg, npoints, csCoef_local, work_local);
            PCHIP(z_seg, rho_seg, npoints, rhoCoef_local, work_local); // 注意：若 PCHIP 要求 
            // 写回全局矩阵
            ssp.cpCoef.block(0, start, 4, npoints) = cpCoef_local;
            ssp.csCoef.block(0, start, 4, npoints) = csCoef_local;
            ssp.rhoCoef.block(0, start, 4, npoints) = rhoCoef_local;
        }

        // Step 3: 如果使用 Cubic Spline 插值
        if (SSPType == SSP_Mode::MODE_S_cCubic)
        {
            int IBCBeg = 0, IBCEnd = 0;

            Eigen::VectorXd z_seg = ssp.z.segment(start, npoints);
            Eigen::VectorXcd cp_val = ssp.cp.segment(start, npoints);
            Eigen::VectorXcd cs_val = ssp.cs.segment(start, npoints);
            Eigen::VectorXcd rho_val = ssp.rho.segment(start, npoints).cast<std::complex<double>>();
            
            // 局部样条系数矩阵（假设 4 行）
            Eigen::MatrixXcd cpSpline_local(4, npoints);
            Eigen::MatrixXcd csSpline_local(4, npoints);
            Eigen::MatrixXcd rhoSpline_local(4, npoints);
            // 第0行初始化为函数值
            cpSpline_local.row(0) = cp_val;
            csSpline_local.row(0) = cs_val;
            rhoSpline_local.row(0) = rho_val;

            // 调用 CSpline（传局部矩阵）
            CSpline(z_seg, cpSpline_local, npoints, IBCBeg, IBCEnd, npoints);
            CSpline(z_seg, csSpline_local, npoints, IBCBeg, IBCEnd, npoints);
            CSpline(z_seg, rhoSpline_local, npoints, IBCBeg, IBCEnd, npoints);

            // 写回全局
            ssp.cpSpline.block(0, start, 4, npoints) = cpSpline_local;
            ssp.csSpline.block(0, start, 4, npoints) = csSpline_local;
            ssp.rhoSpline.block(0, start, 4, npoints) = rhoSpline_local;
            
        }
    }
}

void UpdateHSLoss(double &freq, double &freq0,  Atten_Mode &AttenUnit, HSInfo &HSTop, HSInfo &HSBot)
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
