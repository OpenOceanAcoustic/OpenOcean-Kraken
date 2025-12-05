#include "field.h"

#include <Eigen/Dense>
#include <cmath>
#include <complex>

// 假设输入输出变量类型：
// - M: 模态数量
// - P: 压力场矩阵 (Nz x Nr，复数)
// - C, k: 模态系数和波数向量 (长度M，复数)
// - Option: 计算选项字符串
// - phi: 模态形状矩阵 (M x Nz，复数)
// - ro: 偏移距离向量 (长度Nz，实数)
// - r: 距离向量 (长度Nr，实数)
// - Nz, Nr: 深度和距离网格数量
// - MinExp, TINY: 数值计算参数（根据实际定义补充）

void Evaluate(EigenFunction &eigenfun, EigenParams &eigen, parameters &params,int isz,
     std::complex<float> *uAllSources,
     std::complex<float> *uAllSources_vr,
     std::complex<float> *uAllSources_vz)
{
    // 如果没有模态，返回零压力场
    if (eigen.M <= 0)
    {
        return;
    }

    VectorXcd col_vec;
    col_vec = eigenfun.phiS.col(isz);
    double omega = 2 * pi * params.freqinfo->freq;
    double rho = 1.0; // 假设水的密度为1 g/cm³ （不知道如何导入密度，先在这里设置一个标准值）
    // cout<< "col_vec:\n" << col_vec.real()<<endl;

    // 初始化因子和常数向量
    std::complex<double> factor = I1D * std::sqrt(2.0 * pi) * std::exp(I1D * pi / 4.0);
    Eigen::VectorXcd constants(eigen.M);
    Eigen::VectorXcd constants_vr(eigen.M);
    Eigen::VectorXcd constants_vz(eigen.M);

    // 根据选项计算常数向量
    if (params.SourceType == Source_Mode::MODE_X_Line) // Cylindrical coordinates
    {
        constants = factor * col_vec.array() / eigen.k.array();
        constants_vr = factor * col_vec.array() * eigen.k.array() / (omega * rho);
        constants_vz = factor * col_vec.array() / ( eigen.k.array() * omega * rho * I1D );
    }
    else
    {
        constants = factor * col_vec.array() / eigen.k.array().sqrt();
        constants_vr = factor * col_vec.array() * eigen.k.array().sqrt() / (omega * rho);
        constants_vz = factor / I1D * col_vec.array() / eigen.k.array().sqrt() /(omega * rho);
    }

    // 计算ik向量（波数相关项）
    Eigen::VectorXcd ik = -I1D * eigen.k.array();
    if (params.coherenceType == CoherenceType::Incoherent)
    { // 非相干情况（取实部）
        ik = ik.real().cast<std::complex<double>>();
    }

    // 初始化Cmat矩阵 (M x Nz)：预计算深度相关项
    Eigen::MatrixXcd Cmat(eigen.M, params.Pos->NRz);
    Eigen::MatrixXcd Cmat_vr(eigen.M, params.Pos->NRz);
    Eigen::MatrixXcd Cmat_vz(eigen.M, params.Pos->NRz);
    for (int iz = 0; iz < params.Pos->NRz; ++iz)
    {                                                                 // 0-based索引
        Eigen::VectorXcd exp_terms = ik.array() * params.Pos->Ro(iz); // ik * Rz(iz)
        exp_terms = exp_terms.array().exp();                          // e^(ik * Rz(iz))
        Cmat.col(iz) = constants.array() * eigenfun.phiR.col(iz).array() * exp_terms.array();
        Cmat_vr.col(iz) = constants_vr.array() * eigenfun.phiR.col(iz).array() * exp_terms.array();
        Cmat_vz.col(iz) = constants_vz.array() * eigenfun.dphidz.col(iz).array() * exp_terms.array();
    }

    // 遍历所有距离点计算压力场
    for (int ir = 0; ir < params.Pos->NRr; ++ir)
    { // 距离循环
        // 计算Hankel函数相关项（指数部分）
        Eigen::VectorXcd Hank(eigen.M);
        for (int m = 0; m < eigen.M; ++m)
        { // 模态循环
            // 处理下溢（可选，根据需要启用）
            // if (std::real(ik(m) * params.Pos->Rr(ir)) > MinExp) {
            Hank(m) = std::exp(ik(m) * params.Pos->Rr(ir));
            // } else {
            //     Hank(m) = 0.0;
            // }
        }

        // 计算压力场（相干/非相干情况）
        if (params.coherenceType == CoherenceType::Coherent)
        { // 相干情况
            for (int iz = 0; iz < params.Pos->NRz; ++iz)
            { // 深度循环
                // complex<double> data = (Cmat.col(iz).array() * Hank.array()).sum();
                VectorXcd col_vec1 = Cmat.col(iz);
                complex<double> data = (col_vec1.array() * Hank.array()).sum();
                uAllSources[GetFieldAddr(isz, iz, ir, &params.Pos[0])] = data;

                // 计算水平振速vr
                VectorXcd col_vec_vr = Cmat_vr.col(iz);
                complex<double> data_vr = (col_vec_vr.array() * Hank.array()).sum();
                uAllSources_vr[GetFieldAddr(isz, iz, ir, &params.Pos[0])] = data_vr;

                // 计算垂直振速vz
                VectorXcd col_vec_vz = Cmat_vz.col(iz);
                complex<double> data_vz = (col_vec_vz.array() * Hank.array()).sum();
                uAllSources_vz[GetFieldAddr(isz, iz, ir, &params.Pos[0])] = data_vz;
            }
        }
        else
        { // 非相干情况（取模平方和的平方根）
            for (int iz = 0; iz < params.Pos->NRz; ++iz)
            {
                Eigen::VectorXcd temp = Cmat.col(iz).array() * Hank.array();
                uAllSources[GetFieldAddr(isz, iz, ir, &params.Pos[0])] = temp.array().abs2().sum();                                         // 模平方和
                uAllSources[GetFieldAddr(isz, iz, ir, &params.Pos[0])] = std::sqrt(uAllSources[GetFieldAddr(isz, iz, ir, &params.Pos[0])]); // 开平方

                // 计算水平振速vr
                Eigen::VectorXcd temp_vr = Cmat_vr.col(iz).array() * Hank.array();
                uAllSources_vr[GetFieldAddr(isz, iz, ir, &params.Pos[0])] = temp_vr.array().abs2().sum();                                         // 模平方和
                uAllSources_vr[GetFieldAddr(isz, iz, ir, &params.Pos[0])] = std::sqrt(uAllSources_vr[GetFieldAddr(isz, iz, ir, &params.Pos[0])]); // 开平方

                // 计算垂直振速vz
                Eigen::VectorXcd temp_vz = Cmat_vz.col(iz).array() * Hank.array();
                uAllSources_vz[GetFieldAddr(isz, iz, ir, &params.Pos[0])] = temp_vz.array().abs2().sum();                                         // 模平方和
                uAllSources_vz[GetFieldAddr(isz, iz, ir, &params.Pos[0])] = std::sqrt(uAllSources_vz[GetFieldAddr(isz, iz, ir, &params.Pos[0])]); // 开平方
            }
        }

        // 可选：加入柱面扩展因子
        if (params.SourceType == Source_Mode::MODE_R_Point)
        {
            for (int iz = 0; iz < params.Pos->NRz; ++iz)
            {
                double denom = params.Pos->Rr(ir) + params.Pos->Ro(iz);
                if (std::abs(denom) > 1e-3)
                { // 避免除零
                    uAllSources[GetFieldAddr(isz, iz, ir, &params.Pos[0])] /= std::sqrt(denom);
                    uAllSources_vr[GetFieldAddr(isz, iz, ir, &params.Pos[0])] /= std::sqrt(denom);
                    uAllSources_vz[GetFieldAddr(isz, iz, ir, &params.Pos[0])] /= std::sqrt(denom);
                }
            }
        }
    }
}

void field(EigenFunction &eigenfun, EigenParams &eigen, parameters &params, std::complex<float> *uAllSources, int isz)
{
    VectorXcd constt, sumk;
    constt.resize(eigen.M);
    sumk = VectorXcd::Zero(eigen.M);

    for (int i = 0; i < eigen.M; i++)
    {
        constt(i) = I1D * std::sqrt(2.0 * pi) * std::exp(I1D * pi / 4.0) * eigenfun.phiS(i, isz);
    }

    VectorXcd Hank(eigen.M);
    // VectorXcd phi(eigen.M);

    for (int irr = 0; irr < params.Pos->NRr; irr++)
    {
        // if (irr > 0)
        //     rLeft = std::max(params.Pos->Rr(irr - 1), 0.0);
        // else
        //     rLeft = 0.0;

        double Rr = params.Pos->Rr(irr);
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
        for (int irz = 0; irz < params.Pos->NRz; irz++)
        {
            size_t base = GetFieldAddr(isz, irz, irr, &params.Pos[0]);
            for (int i = 0; i < eigen.M; i++)
            {
                uAllSources[base] += eigenfun.phiR(i, irz) * Hank(i);
            }
        }
    }
}

void export_shd(std::string filename, parameters &params, std::complex<float> *uAllSources)
{
    // std::cout << "开始导出SHD文件: " << filename << ", 数据类型: " << dataType << std::endl;
    int LRecl;
    // 打开文件
    std::ofstream SHDFile = std::ofstream(filename + ".shd", std::ios::binary);
    if (!SHDFile.is_open())
    {
        std::cerr << "无法打开SHD文件: " << filename << ".shd" << std::endl;
        return;
    }
    // std::cout << "SHD文件已打开" << std::endl;

    float Atten = 0.0;

    int Ntheta = 1, NSx = 1, NSy = 1, Nfreq = 1;
    std::vector<double> freqvec(Nfreq, params.freqinfo->freq);
    std::vector<float> theta(Ntheta, 0.0), Sx(NSx, 0.0), Sy(NSy, 0.0);
    std::vector<float> Sz(params.Pos->NSz);
    for (int i = 0; i < params.Pos->NSz; i++)
    {
        Sz[i] = params.Pos->Sz[i];
    }
    std::vector<float> Rz(params.Pos->NRz);
    for (int i = 0; i < params.Pos->NRz; i++)
    {
        Rz[i] = params.Pos->Rz[i];
    }
    std::vector<float> Rr(params.Pos->NRr);
    for (int i = 0; i < params.Pos->NRr; i++)
    {
        Rr[i] = params.Pos->Rr[i];
    }

    std::string plottype;
    switch (params.Pos->GridType)
    {
    case Grid_Mode::MODE_I_Irregular:
        plottype = "irregular ";
        break;
    case Grid_Mode::MODE_R_Rectangular:
        plottype = "rectilin  ";
        break;
    }

    LRecl = std::max(std::max(std::max(std::max(std::max(41, 2 * params.freqinfo->Nfreq), Ntheta), params.Pos->NSz), params.Pos->NRz), 2 * params.Pos->NRr);

    SHDFile.write(reinterpret_cast<char *>(&LRecl), sizeof(int));
    SHDFile.write(params.Title.c_str(), params.Title.size());
    SHDFile.seekp(1 * 4 * LRecl, std::ios::beg);
    SHDFile.write(plottype.c_str(), 80);
    SHDFile.seekp(2 * 4 * LRecl, std::ios::beg);
    SHDFile.write(reinterpret_cast<char *>(&Nfreq), sizeof(int));
    SHDFile.write(reinterpret_cast<char *>(&Ntheta), sizeof(int));
    SHDFile.write(reinterpret_cast<char *>(&NSx), sizeof(int));
    SHDFile.write(reinterpret_cast<char *>(&NSy), sizeof(int));
    SHDFile.write(reinterpret_cast<char *>(&params.Pos->NSz), sizeof(int));
    SHDFile.write(reinterpret_cast<char *>(&params.Pos->NRz), sizeof(int));
    SHDFile.write(reinterpret_cast<char *>(&params.Pos->NRr), sizeof(int));
    float freq = static_cast<float>(params.freqinfo->freq); // 转换double为float
    SHDFile.write(reinterpret_cast<char *>(&freq), sizeof(float));
    SHDFile.write(reinterpret_cast<char *>(&Atten), sizeof(float));
    SHDFile.seekp(3 * 4 * LRecl, std::ios::beg);
    SHDFile.write(reinterpret_cast<char *>(freqvec.data()), Nfreq * sizeof(double));
    SHDFile.seekp(4 * 4 * LRecl, std::ios::beg);
    SHDFile.write(reinterpret_cast<char *>(theta.data()), Ntheta * sizeof(float));
    SHDFile.seekp(5 * 4 * LRecl, std::ios::beg);
    SHDFile.write(reinterpret_cast<char *>(Sx.data()), params.Pos->NSx * sizeof(float));
    SHDFile.seekp(6 * 4 * LRecl, std::ios::beg);
    SHDFile.write(reinterpret_cast<char *>(Sy.data()), params.Pos->NSy * sizeof(float));
    SHDFile.seekp(7 * 4 * LRecl, std::ios::beg);
    SHDFile.write(reinterpret_cast<char *>(Sz.data()), params.Pos->NSz * sizeof(float));
    SHDFile.seekp(8 * 4 * LRecl, std::ios::beg);
    SHDFile.write((char *)Rz.data(), Rz.size() * sizeof(float));
    SHDFile.seekp(9 * 4 * LRecl, std::ios::beg);
    SHDFile.write(reinterpret_cast<char *>(Rr.data()), Rr.size() * sizeof(float));

    // std::cout << "开始写入数据，NSz=" << params.Pos->NSz << ", NRz_per_range=" << params.Pos->NRz_per_range << ", NRr=" << params.Pos->NRr << std::endl;
    // 遍历所有位置并写入声压值的实部和虚部
    for (int isz = 0; isz < params.Pos->NSz; ++isz)
    {
        for (int irz = 0; irz < params.Pos->NRz_per_range; ++irz)
        {
            int recnum = 10 + isz * params.Pos->NRz_per_range + irz;
            SHDFile.seekp(recnum * 4 * LRecl, std::ios::beg);
            for (int ir = 0; ir < params.Pos->NRr; ++ir)
            {
                std::complex<float> &P = uAllSources[GetFieldAddr(isz, irz, ir, params.Pos)];
                SHDFile.write(reinterpret_cast<const char *>(&P), sizeof(P));
            }
        }
        // // 添加进度信息
        // if (isz % 10 == 0 || isz == params.Pos->NSz - 1) {
        //     std::cout << "进度: " << (isz + 1) << "/" << params.Pos->NSz << std::endl;
        // }
    }

    SHDFile.close();
    // std::cout << "SHD文件导出完成: " << filename << std::endl;
}
