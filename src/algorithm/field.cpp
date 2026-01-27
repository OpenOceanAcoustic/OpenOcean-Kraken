#include "field.h"

namespace OpenOceanKraken
{
    void Evaluate(EigenParams &eigen, const OOK_parameters &params, int isz,int iprof,
                  std::complex<float> *uAllSources,
                  std::complex<float> *v_AllSources,
                  std::complex<float> *h_AllSources)
    {
        // 如果没有模态，返回零压力场
        if (eigen.M <= 0)return;


        Eigen::VectorXcd krm = eigen.k.segment(0, eigen.M);
        Eigen::MatrixXcd PsiR = eigen.PsiR.block(0, 0, eigen.M, eigen.PsiR.cols());
        Eigen::MatrixXcd dPsiRdz = eigen.dPsidzR.block(0, 0, eigen.M, eigen.dPsidzR.cols());
        Eigen::MatrixXcd PsiS = eigen.PsiS.block(0, 0, eigen.M, eigen.PsiS.cols());

        Eigen::VectorXcd col_vec;
        col_vec = PsiS.col(isz);
        double c0,rho;
        getSourceEnv(params.SSP.at(iprof),params.Pos.Sz(isz),rho,c0); //计算声源位置的声速和密度
        double omega = 2 * pi * params.freqinfo.freq;


        // 初始化因子和常数向量
        std::complex<double> factor = I1D * std::sqrt(2.0 * pi) * std::exp(I1D * pi / 4.0);
        Eigen::VectorXcd constants(eigen.M);
        Eigen::VectorXcd constants_vr(eigen.M);
        Eigen::VectorXcd constants_vz(eigen.M);

        // 根据选项计算常数向量
        if (params.SourceType == Source_Mode::MODE_X_Line) // Cylindrical coordinates
        {
            constants = factor * col_vec.array() / krm.array();
            constants_vr = factor * col_vec.array() * krm.array() / omega * c0;
            constants_vz = factor * col_vec.array() * c0 / (krm.array() * omega * I1D);
        }
        else
        {
            constants = factor * col_vec.array() / krm.array().sqrt();
            constants_vr = factor * col_vec.array() * krm.array().sqrt() / omega * c0;
            constants_vz = factor * col_vec.array() * c0 / (krm.array().sqrt() * omega * I1D);
        }

        // 计算ik向量（波数相关项）
        Eigen::VectorXcd ik = -I1D * krm.array();
        if (params.coherenceType == CoherenceType::Incoherent)
        { // 非相干情况（取实部）
            ik = ik.real().cast<std::complex<double>>();
        }

        // 初始化Cmat矩阵 (M x Nz)：预计算深度相关项
        Eigen::MatrixXcd Cmat(eigen.M, params.Pos.NRz);
        Eigen::MatrixXcd Cmat_vr(eigen.M, params.Pos.NRz);
        Eigen::MatrixXcd Cmat_vz(eigen.M, params.Pos.NRz);
        for (int iz = 0; iz < params.Pos.NRz; ++iz)
        {                                                                // 0-based索引
            Eigen::VectorXcd exp_terms = ik.array() * params.Pos.Ro(iz); // ik * Ro(iz)
            exp_terms = exp_terms.array().exp();                         // e^(ik * Ro(iz))
            Cmat.col(iz) = constants.array() * PsiR.col(iz).array() * exp_terms.array();
            Cmat_vr.col(iz) = constants_vr.array() * PsiR.col(iz).array() * exp_terms.array();
            Cmat_vz.col(iz) = constants_vz.array() * dPsiRdz.col(iz).array() * exp_terms.array();
        }

        // 遍历所有距离点计算压力场
        for (int ir = 0; ir < params.Pos.NRr; ++ir)
        { // 距离循环
            // 计算Hankel函数相关项（指数部分）
            Eigen::VectorXcd Hank(eigen.M);
            for (int m = 0; m < eigen.M; ++m)
            { // 模态循环
                // 处理下溢（可选，根据需要启用）
                // if (std::real(ik(m) * params.Pos->Rr(ir)) > MinExp) {
                Hank(m) = std::exp(ik(m) * params.Pos.Rr(ir));
                // } else {
                //     Hank(m) = 0.0;
                // }
            }

            // 计算压力场（相干/非相干情况）
            if (params.coherenceType == CoherenceType::Coherent)
            { // 相干情况
                for (int iz = 0; iz < params.Pos.NRz; ++iz)
                { // 深度循环
                    // complex<double> data = (Cmat.col(iz).array() * Hank.array()).sum();
                    Eigen::VectorXcd col_vec1 = Cmat.col(iz);
                    std::complex<double> data = (col_vec1.array() * Hank.array()).sum();
                    uAllSources[GetFieldAddr(isz, iz, ir, &params.Pos)] = data;

                    if (params.is_Velocity)
                    {
                        // 计算水平振速vr
                        Eigen::VectorXcd col_vec_vr = Cmat_vr.col(iz);
                        std::complex<double> data_vr = (col_vec_vr.array() * Hank.array()).sum();
                        v_AllSources[GetFieldAddr(isz, iz, ir, &params.Pos)] = data_vr;

                        // 计算垂直振速vz
                        Eigen::VectorXcd col_vec_vz = Cmat_vz.col(iz);
                        std::complex<double> data_vz = (col_vec_vz.array() * Hank.array()).sum();
                        h_AllSources[GetFieldAddr(isz, iz, ir, &params.Pos)] = data_vz;
                    }
                }
            }
            else
            { // 非相干情况（取模平方和的平方根）
                for (int iz = 0; iz < params.Pos.NRz; ++iz)
                {
                    Eigen::VectorXcd temp = Cmat.col(iz).array() * Hank.array();
                    uAllSources[GetFieldAddr(isz, iz, ir, &params.Pos)] = temp.array().abs2().sum();                                      // 模平方和
                    uAllSources[GetFieldAddr(isz, iz, ir, &params.Pos)] = std::sqrt(uAllSources[GetFieldAddr(isz, iz, ir, &params.Pos)]); // 开平方

                    // 计算水平振速vr
                    if (params.is_Velocity)
                    {
                        Eigen::VectorXcd temp_vr = Cmat_vr.col(iz).array() * Hank.array();
                        v_AllSources[GetFieldAddr(isz, iz, ir, &params.Pos)] = temp_vr.array().abs2().sum();                                    // 模平方和
                        v_AllSources[GetFieldAddr(isz, iz, ir, &params.Pos)] = std::sqrt(v_AllSources[GetFieldAddr(isz, iz, ir, &params.Pos)]); // 开平方

                        // 计算垂直振速vz
                        Eigen::VectorXcd temp_vz = Cmat_vz.col(iz).array() * Hank.array();
                        h_AllSources[GetFieldAddr(isz, iz, ir, &params.Pos)] = temp_vz.array().abs2().sum();                                    // 模平方和
                        h_AllSources[GetFieldAddr(isz, iz, ir, &params.Pos)] = std::sqrt(h_AllSources[GetFieldAddr(isz, iz, ir, &params.Pos)]); // 开平方
                    }
                }
            }

            // 可选：加入柱面扩展因子
            if (params.SourceType == Source_Mode::MODE_R_Point)
            {
                for (int iz = 0; iz < params.Pos.NRz; ++iz)
                {
                    double denom = params.Pos.Rr(ir) + params.Pos.Ro(iz);
                    if (std::abs(denom) > 1e-3)
                    { // 避免除零
                        uAllSources[GetFieldAddr(isz, iz, ir, &params.Pos)] /= std::sqrt(denom);
                        if (params.is_Velocity)
                        {
                            v_AllSources[GetFieldAddr(isz, iz, ir, &params.Pos)] /= std::sqrt(denom);
                            h_AllSources[GetFieldAddr(isz, iz, ir, &params.Pos)] /= std::sqrt(denom);
                        }
                    }
                }
            }
        }
    }

    void field(EigenParams &eigen, const OOK_parameters &params, std::complex<float> *uAllSources, int isz)
    {
        Eigen::VectorXcd constt, sumk;
        constt.resize(eigen.M);
        sumk = Eigen::VectorXcd::Zero(eigen.M);

        for (int i = 0; i < eigen.M; i++)
        {
            constt(i) = I1D * std::sqrt(2.0 * pi) * std::exp(I1D * pi / 4.0) * eigen.PsiS(i, isz);
        }

        Eigen::VectorXcd Hank(eigen.M);
        // VectorXcd phi(eigen.M);

        for (int irr = 0; irr < params.Pos.NRr; irr++)
        {
            // if (irr > 0)
            //     rLeft = std::max(params.Pos->Rr(irr - 1), 0.0);
            // else
            //     rLeft = 0.0;

            double Rr = params.Pos.Rr(irr);
            for (int i = 0; i < eigen.M; i++)
            {
                // sumk(i) += eigen.k(i) * (Rr - rLeft);

                // Hank(i) = constt(i) * std::exp(-I1D * sumk(i));  // coherent   case
                Hank(i) = constt(i) * std::exp(-I1D * eigen.k(i) * Rr); // coherent   case

                if (params.SourceType == Source_Mode::MODE_R_Point) // Cylindrical coordinates
                {
                    if (Rr == 0.0)
                        Hank(i) = 0.0;
                    else
                        Hank(i) /= std::sqrt(eigen.k(i) * Rr);
                }
                else if (params.SourceType == Source_Mode::MODE_X_Line) // Cartesian coordinates
                {
                    Hank(i) = Hank(i) / eigen.k(i);
                }
                else // Scaled cylindrical coordinates
                {
                    Hank(i) = Hank(i) / std::sqrt(eigen.k(i));
                }
            }

            // For each receiver, add up modal contributions
            for (int irz = 0; irz < params.Pos.NRz; irz++)
            {
                size_t base = GetFieldAddr(isz, irz, irr, &params.Pos);
                for (int i = 0; i < eigen.M; i++)
                {
                    uAllSources[base] += eigen.PsiR(i, irz) * Hank(i);
                }
            }
        }
    }

    void getSourceEnv(const ssp::SSPStructure &ssp, const double Sz, double &rho, double &c0)
    {
        // 1. 边界情况：空剖面
        if (ssp.NMedia == 0)
        {
            // 使用顶层半空间（或抛异常，按需）
            rho = ssp.HSTop.rho;
            c0 = ssp.HSTop.alphaR;
            return;
        }

        // 2. 判断 Sz 是否在整体剖面范围内
        double z_min = ssp.z[0];
        double z_max = ssp.z[ssp.z.size() - 1];

        if (Sz <= z_min)
        {
            // 在顶层以上 → 使用 HSTop
            rho = ssp.HSTop.rho;
            c0 = ssp.HSTop.alphaR;
            return;
        }
        if (Sz >= z_max)
        {
            // 在底层以下 → 使用 HSBot
            rho = ssp.HSBot.rho;
            c0 = ssp.HSBot.alphaR;
            return;
        }

        // 3. 定位 Sz 所在的层（线性搜索，因层数通常不多）
        int target_layer = -1;
        for (int i = 0; i < ssp.NMedia; ++i)
        {
            int start = ssp.get_media_start(i);
            int end = ssp.get_media_end(i); // 半开区间 [start, end)
            if (end <= start)
                continue; // 跳过空层

            double layer_z0 = ssp.z[start];
            double layer_z1 = ssp.z[end - 1];

            if (Sz >= layer_z0 && Sz <= layer_z1)
            {
                target_layer = i;
                break;
            }
        }

        // 4. 如果未找到（理论上不应发生，因已检查 z_min/z_max），fallback 到最近点
        if (target_layer == -1)
        {
            // 再次 fallback：用全局最近邻
            Eigen::Index idx;
            (ssp.z.array() - Sz).abs().minCoeff(&idx);
            rho = ssp.rho[idx];
            c0 = ssp.alphaR[idx];
            return;
        }

        // 5. 在目标层内进行一维线性插值（最简单且常用）
        int start = ssp.get_media_start(target_layer);
        int np = ssp.get_media_size(target_layer);

        // 找到 Sz 在该层中的插入位置（lower_bound）
        const double *z_ptr = ssp.z.data() + start;
        auto it = std::lower_bound(z_ptr, z_ptr + np, Sz);
        int local_idx = static_cast<int>(it - z_ptr);

        if (local_idx == 0)
        {
            // 恰好在第一个点
            rho = ssp.rho[start];
            c0 = ssp.alphaR[start];
        }
        else if (local_idx >= np)
        {
            // 恰好在最后一个点（理论上不会，因已检查边界）
            rho = ssp.rho[start + np - 1];
            c0 = ssp.alphaR[start + np - 1];
        }
        else
        {
            // 线性插值 between [local_idx-1, local_idx]
            double z0 = ssp.z[start + local_idx - 1];
            double z1 = ssp.z[start + local_idx];
            double w = (Sz - z0) / (z1 - z0); // 权重 ∈ [0,1]

            rho = ssp.rho[start + local_idx - 1] * (1.0 - w) +
                  ssp.rho[start + local_idx] * w;

            c0 = ssp.alphaR[start + local_idx - 1] * (1.0 - w) +
                 ssp.alphaR[start + local_idx] * w;
        }
    }
}
