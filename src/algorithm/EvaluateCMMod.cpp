// 简正波耦合模块实现文件
#include "EvaluateCMMod.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>

// 构造函数
EvaluateCMMod::EvaluateCMMod() : firstCall(true), iRecProfileR(1) {
}

// 析构函数
EvaluateCMMod::~EvaluateCMMod() {
}

// 主要函数：使用耦合模态理论计算压力场
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
    std::vector<std::complex<double>> A;
    std::string Title;
    
    // 计算使用新剖面的距离（米）
    if (firstCall) {
        for (iProf = nProf; iProf >= 2; iProf--) {
            rProf[iProf] = 500.0 * (rProf[iProf] + rProf[iProf - 1]);
        }
        rProf[nProf + 1] = std::numeric_limits<double>::max();
        firstCall = false;
    }
    
    // 计算模态激励系数 A(mode)
    iProf = 1;
    
    // 获取模态信息
    getModes(fileRoot, iProf, ifreq, MAX_M, rd, nrd, "N", k, phi, M1, freqVec, Nfreq, Title);
    
    // 确保模态数量不超过计算的模态数量
    m = std::min(m, M1);
    
    // 初始化A向量
    A.resize(m);
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
        } else { // 没有新段，只根据距离步长推进相位
            if (ir == 0) { // 第一个距离
                for (int mode = 0; mode < m; mode++) {
                    A[mode] *= std::exp(-I1D * k[mode] * r[ir]);
                }
            } else {
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

// 计算界面左侧的压力场
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

// 计算半空间尾部贡献
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

// 计算Pekeris根
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

// 从文件读取模态数据的辅助函数（需要根据实际文件格式实现）
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