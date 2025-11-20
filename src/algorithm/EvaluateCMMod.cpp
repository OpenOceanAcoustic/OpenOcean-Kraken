//===================================================================
// 简正波耦合模块实现文件 EvaluateCMMod.cpp
// 
// 功能: 实现海洋声学中的简正波耦合理论算法，用于计算水平不均匀海洋环境中的声场
// 实现细节: 
//        - 基于Fortran版本的EvaluateCMMod模块转换为C++实现
//        - 使用C++标准库和项目全局参数定义
//        - 实现简正波在水平不均匀环境中的传播和耦合计算
// 
// 版本: 1.0
// 作者: OpenOcean-Kraken开发团队
//===================================================================
#include "EvaluateCMMod.h"
#include <cmath>          // 用于数学函数运算
#include <fstream>        // 用于文件操作
#include <iostream>       // 用于输入输出
#include <stdexcept>      // 用于异常处理

/**
 * @brief 构造函数
 * 
 * 初始化EvaluateCMMod类的实例，设置所有必要的初始状态变量。
 * 这包括第一次调用标志和记录位置指针，用于后续模态计算和文件读取。
 */
EvaluateCMMod::EvaluateCMMod() : firstCall(true), iRecProfileR(1) {
    // 构造函数初始化工作
}

/**
 * @brief 析构函数
 * 
 * 清理EvaluateCMMod类实例使用的资源。
 * 当前实现中没有动态分配的资源需要释放，但保留此函数以备将来扩展。
 */
EvaluateCMMod::~EvaluateCMMod() {
    // 析构函数清理工作
}

// 主要函数：使用耦合模态理论计算压力场
/**
 * @brief 核心算法实现：使用简正波耦合理论计算声场
 * 
 * 本函数实现了简正波耦合方法的核心算法，通过以下步骤计算水平不均匀海洋环境中的声场：
 * 1. 初始化并处理剖面对应的距离数组
 * 2. 计算初始模态激励系数
 * 3. 沿水平距离传播，处理界面处的模态耦合
 * 4. 在每个接收点计算总声场
 * 
 * 算法原理：
 * - 当声波在水平不均匀海洋环境中传播时，不同位置的简正波模态会发生耦合
 * - 本算法将水平空间划分为多个均匀段，在每个段内声波以简正波形式传播
 * - 在段交界处计算模态间的耦合系数，实现能量从一个模态到另一个模态的转移
 * - 使用耦合系数更新模态振幅，继续向前传播
 */
void EvaluateCMMod::evaluateCM(const std::string& fileRoot, std::vector<double>& rProf, const int& nProf, 
                              const std::vector<std::complex<double>>& phiS, 
                              std::vector<std::vector<std::complex<double>>>& phi, 
                              const std::vector<double>& rd, const int& nrd, 
                              const std::vector<double>& r, const int& nr, 
                              std::vector<std::complex<double>>& k, int& m, 
                              const std::string& option, 
                              std::vector<std::vector<std::complex<double>>>& p) {
    
    // 初始化变量
    int ir, iProf, ird, ifreq = 1, Nfreq, M1;
    std::vector<double> freqVec(MAX_NFREQ);
    std::vector<std::complex<double>> A; // 模态振幅系数数组
    std::string Title; // 模态文件标题信息
    
    // 计算使用新剖面的距离（米）
    if (firstCall) {
        for (iProf = nProf; iProf >= 2; iProf--) {
            rProf[iProf] = 500.0 * (rProf[iProf] + rProf[iProf - 1]);
        }
        rProf[nProf + 1] = std::numeric_limits<double>::max();
        firstCall = false;
    }
    
    // 计算模态激励系数 A(mode)
    iProf = 1;  // 从第一个剖面开始计算
    
    // 获取模态信息
    // 读取第一个剖面对应的模态信息，包括波数和模态函数
    getModes(fileRoot, iProf, ifreq, MAX_M, rd, nrd, "N", k, phi, M1, freqVec, Nfreq, Title);
    
    // 确保模态数量不超过计算的模态数量
    m = std::min(m, M1);
    
    // 初始化A向量（模态振幅系数）
    A.resize(m);  // A[i] 表示第i个模态的激励振幅系数
    
    // 计算相位因子，用于模态激励系数计算
    // 这是简正波理论中的标准相位因子，与格林函数的辐射条件相关
    const std::complex<double> phaseFactor = I1D * pi / 4.0;
    
    // 根据坐标类型计算激励系数
        if (option.substr(0, 1) == "X") { // Cartesian coordinates
            for (int mode = 0; mode < m; mode++) {
                A[mode] = std::sqrt(2.0 * pi) * std::exp(phaseFactor) * phiS[mode] / k[mode];
            }
        } else { // Cylindrical coordinates
            for (int mode = 0; mode < m; mode++) {
                A[mode] = I1D * std::sqrt(2.0 * pi) * std::exp(phaseFactor) * phiS[mode] / std::sqrt(k[mode]);
            }
        }
    
    // 沿距离前进
        for (ir = 0; ir < nr; ir++) {
            if (r[ir] > rProf[iProf + 1]) { // 进入新的距离段?
                iProf++;
                
                // 前进到界面
                if (ir == 0) { // 第一个距离
                    for (int mode = 0; mode < m; mode++) {
                        A[mode] *= std::exp(-I1D * k[mode] * rProf[iProf]);
                    }
                } else {
                    for (int mode = 0; mode < m; mode++) {
                        A[mode] *= std::exp(-I1D * k[mode] * (rProf[iProf] - r[ir - 1]));
                    }
                }
            
            // 这里是跨界面的地方
            if (iProf <= nProf) {
                newProfile(fileRoot, ifreq, k, phi, m, rd, nrd, A);
                std::cout << "New profile read at range " << r[ir] / 1000.0 << ", iProf = " << iProf 
                          << ", #modes = " << m << std::endl;
            }
            
            // 是否有其他段需要跨越? 推进每个段的相位
            while (r[ir] > rProf[iProf + 1]) {
                iProf++;
                for (int mode = 0; mode < m; mode++) {
                    A[mode] *= std::exp(-I1D * k[mode] * (rProf[iProf] - rProf[iProf - 1]));
                }
                
                if (iProf <= nProf) {
                    newProfile(fileRoot, ifreq, k, phi, m, rd, nrd, A);
                }
            }
            
            // 前进经过最后一个界面的剩余距离
            for (int mode = 0; mode < m; mode++) {
                A[mode] *= std::exp(-I1D * k[mode] * (r[ir] - rProf[iProf]));
            }
        } else { // 没有进入新段，只在当前剖面内传播
            if (ir == 0) { // 第一个距离点
                // 从源点到计算点的相位因子
                for (int mode = 0; mode < m; mode++) {
                    A[mode] *= std::exp(-I1D * k[mode] * r[ir]);
                }
            } else { // 后续距离点
                // 从上一个计算点到当前点的相位因子
                for (int mode = 0; mode < m; mode++) {
                    A[mode] *= std::exp(-I1D * k[mode] * (r[ir] - r[ir - 1]));
                }
            }
        }
        
        // 对每个接收器累加模态贡献
        for (ird = 0; ird < nrd; ird++) {
            std::complex<double> sum(0.0, 0.0);
            for (int mode = 0; mode < m; mode++) {
                sum += A[mode] * phi[mode][ird];
            }
            
            if (option.substr(0, 1) == "R" && r[ir] != 0.0) {
                p[ird][ir] = sum / std::sqrt(r[ir]);
            } else {
                p[ird][ir] = sum;
            }
        }
    }
}

// 处理新剖面对应的模态
void EvaluateCMMod::newProfile(const std::string& fileRoot, const int& ifreq, 
                              std::vector<std::complex<double>>& k, 
                              std::vector<std::vector<std::complex<double>>>& phiR, 
                              int& mr, const std::vector<double>& rd, const int& nrd, 
                              std::vector<std::complex<double>>& a) {
    
    // 变量初始化
    int mode, nr, nTot, ml, ir, iz;
    std::vector<int> irdVec(nrd);
    double depthTL, depthBL, depthTR, depthBR, rhoBR, rhoTR, zt;
    std::vector<std::complex<double>> p, phi, phiTmp;
    std::vector<std::complex<double>> gamTL(mr), gamBL(mr), phiTL(mr), phiBL(mr);
    std::complex<double> gamTR, gamBR, phiTR(0.0), phiBR(0.0), kTop2R, kBot2R;
    std::complex<double> gamma2;
    char bcTopR, bcBotR;
    std::vector<double> z, w(nrd);
    
    // 计算界面左侧的压力
    pLeft(fileRoot, ifreq, iRecProfileR, a, k, z, mr, p, nr, nTot, 
          bcTopR, rhoTR, kTop2R, depthTR, bcBotR, rhoBR, kBot2R, depthBR, 
          gamTL, gamBL, depthTL, depthBL, phiTL, phiBL, ml);
    
    // 计算接收器深度处模态插值的权重
    weight(z, nTot, rd, nrd, w, irdVec);
    
    // 分配内存
    phi.resize(nr);
    phiTmp.resize(nTot);
    
    // *** 将压力 P 投影到新段中的每个模态 ***
    for (mode = 0; mode < mr; mode++) {
        // 这里需要从文件读取模态，暂时使用占位符
        // READ(ModeFile, REC = IRecProfile + 1 + mode) (phi(iz), iz = 1, Nr)
        
        // 模拟模态读取，实际实现需要从文件读取
        for (iz = 0; iz < nr; iz++) {
            // 实际应用中需要替换为真实的模态值
            phi[iz] = std::complex<double>(0.0, 0.0);
        }
        
        if (bcTopR == 'A') {
            phiTR = phi[0];
            gamma2 = k[mode] * k[mode] - kTop2R;
            gamTR = pekerisRoot(gamma2);
        }
        
        if (bcBotR == 'A') {
            phiBR = phi[nr - 1];
            gamma2 = k[mode] * k[mode] - kBot2R;
            gamBR = pekerisRoot(gamma2);
        }
        
        // 在先前段的网格上制表新模态
        for (iz = 0; iz < nTot; iz++) {
            zt = z[iz];
            if (zt > depthBR) {
                if (bcBotR == 'A') {
                    phiTmp[iz] = phiBR * std::exp(-gamBR * (zt - depthBR));
                }
            } else if (zt < depthTR) {
                if (bcTopR == 'A') {
                    phiTmp[iz] = phiTR * std::exp(-gamTR * (depthTR - zt));
                }
            } else {
                phiTmp[iz] = phi[iz];
            }
        }
        
        // 计算新振幅: A = Integral[ P(z) * phi(z) dz ]
        std::complex<double> sum1(0.0, 0.0);
        for (iz = 0; iz < nTot; iz++) {
            sum1 += p[iz] * phiTmp[iz];
        }
        
        // 上半空间的贡献
        if (bcTopR == 'A') {
            sum1 += calculateTail(z[0], phiTL, gamTL, depthTL, ml, phiTR / rhoTR, gamTR, depthTR);
        }
        
        // 下半空间的贡献
        if (bcBotR == 'A') {
            sum1 += calculateTail(z[nTot - 1], phiBL, gamBL, depthBL, ml, phiBR / rhoBR, gamBR, depthBR);
        }
        
        a[mode] = sum1;
        
        // 在接收器深度处制表模态
        for (ir = 0; ir < nrd; ir++) {
            phiR[mode][ir] = std::complex<double>(0.0, 0.0);
            if (rd[ir] < depthTR) { // 接收器在上半空间
                if (bcTopR == 'A') {
                    phiR[mode][ir] = phiTR * std::exp(-gamTR * (depthTR - rd[ir]));
                }
            } else if (rd[ir] > depthBR) { // 接收器在下半空间
                if (bcBotR == 'A') {
                    phiR[mode][ir] = phiBR * std::exp(-gamBR * (rd[ir] - depthBR));
                }
            } else { // 接收器在内部区域
                iz = irdVec[ir];
                phiR[mode][ir] = phi[iz] + w[ir] * (phi[iz + 1] - phi[iz]);
            }
        }
    }
}

/**
 * @brief 计算界面左侧的压力场
 * 
 * 从文件中读取前一个剖面对应的模态数据，计算界面左侧的压力场分布。
 * 这是模态耦合计算中的关键步骤，用于获取界面处的声场信息。
 * 
 * @param fileRoot 模态文件路径前缀
 * @param ifreq 频率索引
 * @param iRecProfile 记录位置指针（输入/输出）
 * @param a 模态振幅系数（输入/输出）
 * @param k 波数向量（输出）
 * @param z 深度坐标数组（输出）
 * @param m 模态数量（输出）
 * @param p 压力场数组（输出）
 * @param nr 深度点数（输出）
 * @param nTot 总深度点数（输出）
 * @param bcTop 顶部边界条件（输出）
 * @param rhoT 顶部介质密度（输出）
 * @param kTop2 顶部介质波数平方（输出）
 * @param depthT 顶部深度（输出）
 * @param bcBot 底部边界条件（输出）
 * @param rhoB 底部介质密度（输出）
 * @param kBot2 底部介质波数平方（输出）
 * @param depthB 底部深度（输出）
 * @param gamTL 顶部半空间衰减常数（输出）
 * @param gamBL 底部半空间衰减常数（输出）
 * @param depthTL 左侧顶部深度（输出）
 * @param depthBL 左侧底部深度（输出）
 * @param phiTL 左侧顶部模态值（输出）
 * @param phiBL 左侧底部模态值（输出）
 * @param ml 左侧模态数量（输出）
 */
void EvaluateCMMod::pLeft(const std::string& fileRoot, const int& ifreq, int& iRecProfile, 
                         std::vector<std::complex<double>>& a, std::vector<std::complex<double>>& k, 
                         std::vector<double>& z, int& m, std::vector<std::complex<double>>& p, 
                         int& nr, int& nTot, char& bcTop, double& rhoT, std::complex<double>& kTop2, 
                         double& depthT, char& bcBot, double& rhoB, std::complex<double>& kBot2, 
                         double& depthB, std::vector<std::complex<double>>& gamTL, 
                         std::vector<std::complex<double>>& gamBL, double& depthTL, double& depthBL, 
                         std::vector<std::complex<double>>& phiTL, std::vector<std::complex<double>>& phiBL, 
                         int& ml) {
    
    // 初始化变量
    int nMat, med, mode, nMedia, iz, izL, nfreq, lRecl, mr;
    std::vector<int> n, nL;
    std::vector<double> zL, depthL, depth, rhoL, rho, freqVec(MAX_NFREQ);
    double dbBelow, rhomed, rhoBel, h;
    std::complex<double> cPT, cPB, cST, cSB, gamma2;
    std::vector<std::complex<double>> phi, pL;
    std::string title;
    std::vector<std::string> material, materialL;
    // 使用已定义的I1D常量作为虚数单位
    
    // 读取上一段末尾的模态信息
    iRecProfile = iRecProfileR;
    
    // 读取模态头部信息
    readModeHeader(fileRoot, iProf - 1, iRecProfile, lRecl, title, freqVec, nfreq, nMedia, 
                  nL[0], nMat, n, materialL, depthL, rhoL, zL);
    
    // 计算界面处的压力
    phi.resize(nL[0]);
    pL.resize(nL[0], std::complex<double>(0.0, 0.0));
    
    // 读取特征值 k(I)
    readWavenumbers(iRecProfile, ifreq, k, ml, MAX_M, lRecl);
    if (ml == 0) return;
    
    // 读取顶部和底部半空间信息
    // 实际应用中需要从文件读取
    bcTop = 'A'; // 假设值，实际需要从文件读取
    bcBot = 'A'; // 假设值，实际需要从文件读取
    cPT = std::complex<double>(1500.0, 0.0); // 假设值
    cPB = std::complex<double>(1600.0, 0.0); // 假设值
    rhoT = 1000.0; // 假设值
    rhoB = 1800.0; // 假设值
    depthTL = 0.0; // 假设值
    depthBL = 1000.0; // 假设值
    
    // 设置 kTop2, kBot2 基于频率
    if (bcTop == 'A') {
        double omega = 2.0 * PI * freqVec[ifreq];
        kTop2 = (omega / cPT) * (omega / cPT);
    }
    if (bcBot == 'A') {
        double omega = 2.0 * PI * freqVec[ifreq];
        kBot2 = (omega / cPB) * (omega / cPB);
    }
    
    // 模态求和计算左侧的压力 pL
    int modeLimit = std::min(ml, m);
    gamTL.resize(modeLimit);
    gamBL.resize(modeLimit);
    phiTL.resize(modeLimit);
    phiBL.resize(modeLimit);
    
    for (mode = 0; mode < modeLimit; mode++) {
        // 从文件读取模态，暂时使用占位符
        // READ(ModeFile, REC = IRecProfile + 1 + mode) phi
        for (iz = 0; iz < nL[0]; iz++) {
            // 实际应用中需要替换为真实的模态值
            phi[iz] = std::complex<double>(0.0, 0.0);
        }
        
        // 累加模态贡献
        for (iz = 0; iz < nL[0]; iz++) {
            pL[iz] += a[mode] * phi[iz];
        }
        
        // 半空间信息
        phiTL[mode] = a[mode] * phi[0];
        phiBL[mode] = a[mode] * phi[nL[0] - 1];
        gamTL[mode] = std::complex<double>(0.0, 0.0);
        gamBL[mode] = std::complex<double>(0.0, 0.0);
        
        if (bcTop == 'A') { // 顶部半空间
            gamma2 = k[mode] * k[mode] - kTop2;
            gamTL[mode] = pekerisRoot(gamma2);
        }
        
        if (bcBot == 'A') { // 底部半空间
            gamma2 = k[mode] * k[mode] - kBot2;
            gamBL[mode] = pekerisRoot(gamma2);
        }
    }
    
    // 设置记录指针到下一组模态的开始
    iRecProfileR = iRecProfile + ml + 3 + (2 * ml - 1) / lRecl;
    
    // 读取新段中的模态数据
    iRecProfile = iRecProfileR;
    
    // 读取模态头部信息
    readModeHeader(fileRoot, iProf, iRecProfile, lRecl, title, freqVec, nfreq, nMedia, 
                  nr, nMat, n, material, depth, rho, z);
    
    // 读取特征值 k(I)
    readWavenumbers(iRecProfile, ifreq, k, mr, MAX_M, lRecl);
    if (mr == 0) return;
    m = mr;
    
    // 读取顶部和底部半空间信息
    // 实际应用中需要从文件读取
    bcTop = 'A'; // 假设值
    bcBot = 'A'; // 假设值
    cPT = std::complex<double>(1500.0, 0.0); // 假设值
    cPB = std::complex<double>(1600.0, 0.0); // 假设值
    rhoT = 1000.0; // 假设值
    rhoB = 1800.0; // 假设值
    depthT = 0.0; // 假设值
    depthB = 1000.0; // 假设值
    
    if (z[0] != depthT || z[nr - 1] != depthB) {
        std::cerr << "Fatal Error: modes must be tabulated throughout the ocean and sediment to compute the coupling coefs." << std::endl;
        std::cerr << "depths: " << depthT << ", " << depthB << std::endl;
        std::cerr << "z: " << z[0] << ", " << z[nr - 1] << std::endl;
        throw std::runtime_error("Incompatible depth values");
    }
    
    // 上倾? 用zL的数据扩展z向量
    nTot = nr;
    for (izL = 0; izL < nL[0]; izL++) {
        if (zL[izL] > z[nTot - 1]) {
            z.push_back(zL[izL]);
            nTot++;
        }
    }
    
    // 在新网格上重新制表压力
    izL = 0;
    med = 0;
    rhomed = rho[0];
    
    // 此界面下方的下一个界面深度
    if (med < nMedia) {
        dbBelow = depth[med + 1];
        rhoBel = rho[med + 1];
    } else {
        dbBelow = depthB;
        rhoBel = rhoB;
    }
    
    // 分配p的内存
    p.resize(nTot, std::complex<double>(0.0, 0.0));
    
    // 网格点循环
    for (iz = 0; iz < nTot; iz++) {
        double zt = z[iz];
        
        // 获取介质密度
        if (zt < depthT) {
            rhomed = rhoT;
        } else if (zt > depthB) {
            rhomed = rhoB;
        } else if (med < nMedia) {
            if (zt > dbBelow) {
                med++;
                rhomed = rho[med];
            }
        }
        
        // 此界面下方的下一个界面深度
        if (med < nMedia) {
            dbBelow = depth[med + 1];
            rhoBel = rho[med + 1];
        } else {
            dbBelow = depthB;
            rhoBel = rhoB;
        }
        
        while (zt > zL[izL + 1] && izL < nL[0] - 1) {
            izL++;
        }
        
        // 计算该深度的P
        if (zt > depthBL) { // 下半空间
            if (bcBot == 'A') {
                std::complex<double> sum(0.0, 0.0);
                for (mode = 0; mode < ml; mode++) {
                    sum += phiBL[mode] * std::exp(-gamBL[mode] * (zt - depthBL));
                }
                p[iz] = sum;
            }
        } else if (zt < depthTL) { // 上半空间
            if (bcTop == 'A') {
                std::complex<double> sum(0.0, 0.0);
                for (mode = 0; mode < ml; mode++) {
                    sum += phiTL[mode] * std::exp(-gamTL[mode] * (depthTL - zt));
                }
                p[iz] = sum;
            }
        } else { // 线性插值
            p[iz] = pL[izL] + 
                   (zt - zL[izL]) / (zL[izL + 1] - zL[izL]) * 
                   (pL[izL + 1] - pL[izL]);
        }
        
        // 计算网格宽度 h
        if (iz == 0) { // 第一个点
            h = 0.5 * (z[1] - z[0]) / rhomed;
        } else if (iz == nTot - 1) { // 最后一个点
            h = 0.5 * (z[nTot - 1] - z[nTot - 2]) / rhomed;
        } else { // 恰好在界面上方或下方的点
            if (z[iz - 1] < dbBelow && z[iz + 1] > dbBelow) {
                h = 0.5 * (z[iz + 1] / rhoBel - z[iz - 1] / rhomed -
                         dbBelow / rhoBel + dbBelow / rhomed);
            } else {
                h = 0.5 * (z[iz + 1] - z[iz - 1]) / rhomed;
            }
        }
        
        p[iz] *= h;
    }
}

/**
 * @brief 计算半空间尾部贡献
 * 
 * 计算半空间中声场的尾部贡献，用于处理声场在无限半空间中的传播特性。
 * 该函数考虑了界面两侧半空间中声场的耦合效应。
 * 
 * @param depth 计算点深度
 * @param phiL 左侧半空间模态值
 * @param gamL 左侧半空间衰减常数
 * @param depthL 左侧半空间界面深度
 * @param ml 左侧模态数量
 * @param phiR 右侧半空间模态值
 * @param gamR 右侧半空间衰减常数
 * @param depthR 右侧半空间界面深度
 * @return 半空间尾部贡献值
 */
std::complex<double> EvaluateCMMod::calculateTail(const double& depth, const std::vector<std::complex<double>>& phiL, 
                                                const std::vector<std::complex<double>>& gamL, const double& depthL, 
                                                const int& ml, const std::complex<double>& phiR, 
                                                const std::complex<double>& gamR, const double& depthR) {
    
    std::complex<double> fr = phiR * std::exp(-gamR * (depth - depthR));
    std::complex<double> tail(0.0, 0.0);
    
    if (depth == depthL) {
        for (int mode = 0; mode < ml; mode++) {
            tail += phiL[mode] / (gamL[mode] + gamR);
        }
    } else {
        for (int mode = 0; mode < ml; mode++) {
            tail += phiL[mode] * std::exp(-gamL[mode] * (depth - depthL)) / (gamL[mode] + gamR);
        }
    }
    
    return fr * tail;
}

/**
 * @brief 计算Pekeris根
 * 
 * 计算Pekeris波导模型中的波数根，用于处理半空间中的声波传播。
 * 根据gamma2的值计算正确的平方根，确保返回的根具有物理意义（正实部）。
 * 
 * @param gamma2 gamma值的平方，由k^2 - k_0^2计算得到
 * @return 具有物理意义的根值
 */
std::complex<double> EvaluateCMMod::pekerisRoot(const std::complex<double>& gamma2) {
    // 实现Pekeris根计算，这是一个特定于海洋声学的函数
    // 这里提供一个简单的实现，实际应用中可能需要更复杂的处理
    double re = std::real(gamma2);
    double im = std::imag(gamma2);
    
    if (re >= 0 && im == 0) {
        // 实数情况，返回正根
        return std::complex<double>(std::sqrt(re), 0.0);
    } else {
        // 复数情况，返回具有正实部的根
        std::complex<double> root = std::sqrt(gamma2);
        if (std::real(root) < 0) {
            root = -root;
        }
        return root;
    }
}

/**
 * @brief 从文件读取模态数据
 * 
 * 读取指定剖面对应的模态数据，包括波数和模态函数。
 * 这些数据用于声场计算和模态耦合分析。
 * 
 * @param fileRoot 模态文件路径前缀
 * @param iProf 剖面索引
 * @param ifreq 频率索引
 * @param maxM 最大模态数量
 * @param rd 接收器深度数组
 * @param nrd 接收器深度点数
 * @param option 计算选项
 * @param k 波数向量（输出）
 * @param phi 模态函数矩阵（输出）
 * @param m1 实际模态数量（输出）
 * @param freqVec 频率数组（输出）
 * @param nfreq 频率点数（输出）
 * @param title 文件标题（输出）
 */
void EvaluateCMMod::getModes(const std::string& fileRoot, const int& iProf, const int& ifreq, const int& maxM, 
                            const std::vector<double>& rd, const int& nrd, const std::string& option, 
                            std::vector<std::complex<double>>& k, std::vector<std::vector<std::complex<double>>>& phi, 
                            int& m1, std::vector<double>& freqVec, int& nfreq, std::string& title) {
    // 实际应用中需要实现文件读取逻辑
    // 这里提供一个简单的占位实现
    m1 = 10; // 假设10个模态
    k.resize(m1, std::complex<double>(1.0, 0.0));
    phi.resize(m1, std::vector<std::complex<double>>(nrd, std::complex<double>(0.0, 0.0)));
    freqVec.resize(1, 1000.0); // 假设1kHz频率
    nfreq = 1;
    title = "EvaluateCMMod - Mode Data";
}

/**
 * @brief 读取模态文件头部信息
 * 
 * 读取模态文件的头部信息，包括频率、介质参数、深度网格等。
 * 这些信息对于正确解析模态数据和进行声场计算至关重要。
 * 
 * @param fileRoot 模态文件路径前缀
 * @param iProf 剖面索引
 * @param iRecProfile 记录位置指针（输入/输出）
 * @param lRecl 记录长度（输出）
 * @param title 文件标题（输出）
 * @param freqVec 频率数组（输出）
 * @param nfreq 频率点数（输出）
 * @param nMedia 介质数量（输出）
 * @param nl 深度点数（输出）
 * @param nMat 材料数量（输出）
 * @param n 每层介质的网格点数（输出）
 * @param material 材料名称（输出）
 * @param depth 界面深度（输出）
 * @param rho 介质密度（输出）
 * @param z 深度坐标数组（输出）
 */
void EvaluateCMMod::readModeHeader(const std::string& fileRoot, const int& iProf, int& iRecProfile, 
                                  int& lRecl, std::string& title, std::vector<double>& freqVec, 
                                  int& nfreq, int& nMedia, int& nl, int& nMat, 
                                  std::vector<int>& n, std::vector<std::string>& material, 
                                  std::vector<double>& depth, std::vector<double>& rho, 
                                  std::vector<double>& z) {
    // 实际应用中需要实现文件读取逻辑
    // 这里提供一个简单的占位实现
    lRecl = 100; // 假设记录长度
    title = "EvaluateCMMod - Mode Header";
    freqVec.resize(1, 1000.0);
    nfreq = 1;
    nMedia = 2; // 假设2个介质
    nl = 100; // 假设100个深度点
    nMat = 2;
    n.resize(nMat, 50);
    material.resize(nMat, "water");
    depth.resize(nMedia, 500.0);
    rho.resize(nMedia, 1000.0);
    z.resize(nl);
    for (int i = 0; i < nl; i++) {
        z[i] = i * 10.0; // 假设10米间隔
    }
}

/**
 * @brief 读取波数数据
 * 
 * 从模态文件中读取指定频率下的波数数据。
 * 波数是简正波分析中的关键参数，用于计算声场的相位和衰减。
 * 
 * @param iRecProfile 记录位置指针
 * @param ifreq 频率索引
 * @param k 波数向量（输出）
 * @param ml 实际模态数量（输出）
 * @param maxM 最大模态数量
 * @param lRecl 记录长度
 */
void EvaluateCMMod::readWavenumbers(const int& iRecProfile, const int& ifreq, 
                                  std::vector<std::complex<double>>& k, int& ml, const int& maxM, 
                                  const int& lRecl) {
    // 实际应用中需要实现文件读取逻辑
    // 这里提供一个简单的占位实现
    ml = 10; // 假设10个模态
    k.resize(ml);
    for (int i = 0; i < ml; i++) {
        k[i] = std::complex<double>(2.0 + i * 0.1, 0.001);
    }
}

/**
 * @brief 计算接收器深度处模态插值的权重
 * 
 * 计算接收器所在深度处需要的模态函数插值权重。
 * 由于接收器可能位于离散网格点之间，需要通过线性插值获取该深度处的模态值。
 * 
 * @param z 深度网格数组
 * @param nTot 总深度点数
 * @param rd 接收器深度数组
 * @param nrd 接收器数量
 * @param w 插值权重数组（输出）
 * @param ird 接收器所在区间的索引数组（输出）
 */
void EvaluateCMMod::weight(const std::vector<double>& z, const int& nTot, const std::vector<double>& rd, 
                          const int& nrd, std::vector<double>& w, std::vector<int>& ird) {
    // 计算接收器深度处模态插值的权重
    int iz;
    
    for (int ir = 0; ir < nrd; ir++) {
        // 找到深度rd[ir]对应的z中的位置
        for (iz = 0; iz < nTot - 1; iz++) {
            if (rd[ir] >= z[iz] && rd[ir] <= z[iz + 1]) {
                break;
            }
        }
        
        // 确保在有效范围内
        if (iz >= nTot - 1) {
            iz = nTot - 2;
        }
        
        ird[ir] = iz;
        
        // 计算线性插值权重
        if (z[iz + 1] > z[iz]) {
            w[ir] = (rd[ir] - z[iz]) / (z[iz + 1] - z[iz]);
        } else {
            w[ir] = 0.0;
        }
    }
}