#include "init.h"

// 初始化有限差分方程
void Initialize(parameters& params, KrakenMatrix& kramtrx, int ntimes) {
    bool ElasticFlag = false;
    int NPoints = 0;

    double Two_h;
    double cp2, cs2;
    double omega2 = SQ(2 * pi * params.freqinfo->freq);
    
    // 初始化变量
    double& Clow = params.Clow;
    double cMin = 1e8;
    double& cHigh = params.Chigh;
    kramtrx.FirstAcoustic = -1;
    kramtrx.Loc.resize(params.NMedia+1);
    kramtrx.Loc[0] = 0; // C++使用0-based索引
    
    // 计算总网格点数
    kramtrx.N.resize(params.NMedia);
    kramtrx.h.resize(params.NMedia);
    for (int i = 0; i < params.NMedia; ++i) {
        kramtrx.N(i) = params.SSP[i].N*ntimes;
        kramtrx.h(i) = params.SSP[i].depth / kramtrx.N(i);
        NPoints += kramtrx.N(i);
    }
    NPoints += params.NMedia;
    
    // 分配内存（在C++中使用vector会自动管理内存）
    kramtrx.B1.resize(NPoints);
    kramtrx.B1C.resize(NPoints);
    kramtrx.B2.resize(NPoints);
    kramtrx.B3.resize(NPoints);
    kramtrx.B4.resize(NPoints);
    kramtrx.rho.resize(NPoints);
    kramtrx.h.resize(params.NMedia);
    
    // 处理每个介质层
    for (int im = 0; im < params.NMedia; ++im) { // C++使用0-based索引
        // 计算当前层的起始位置
        if (im != 0) {
            kramtrx.Loc(im) = kramtrx.Loc(im - 1) + kramtrx.N(im - 1) + 1;
        }
        
        int ii = kramtrx.Loc(im); // C++使用0-based索引，不需要+1
        
        // 调用EvaluateSSP函数
        SSPStructure SSP = params.SSP[im];
        SSP.N = kramtrx.N(im);

        EvaluateSSP(SSP, params.SSPType);
        
        // 加载有限差分方程的对角线
        if (std::real(SSP.cs[0]) == 0.0) { // 声学介质情况
            SSP.Material = Media_Mode::MODE_A_Acoustic;
            if (kramtrx.FirstAcoustic == -1) {
                kramtrx.FirstAcoustic = im;
            }
            kramtrx.LastAcoustic = im;
            
            // 计算当前层的最小声速
            double min_cp = 1e8;
            for (int j = 0; j < SSP.N + 1; ++j) {
                min_cp = std::min(min_cp, std::real(SSP.cp_int(j)));
            }
            cMin = std::min(cMin, min_cp);
            
            // 计算B1和B1C
            double h2 = SQ(kramtrx.h(im));
            for (int j = 0; j < SSP.N + 1; ++j) {
                cp2 = real(SQ(SSP.cp_int(j)));
                double val = omega2 / cp2;
                kramtrx.B1(ii+j) = -2.0 + h2 * std::real(val);
                kramtrx.B1C(ii+j) = std::imag(val);
            }
        } else { // 弹性介质情况
            SSP.Material = Media_Mode::MODE_E_Elastic;
            ElasticFlag = true;
            Two_h = 2.0 * kramtrx.h[im];
            
            for (int j = 0; j < SSP.N + 1; ++j) {
                cMin = std::min(std::real(SSP.cs_int[j]), cMin);
                
                cp2 = SQ(std::real(SSP.cp_int(j)));
                cs2 = SQ(std::real(SSP.cs_int(j)));
                
                kramtrx.B1(ii+j) = Two_h / (SSP.rho_int(j) * cs2);
                kramtrx.B2(ii+j) = Two_h / (SSP.rho_int(j) * cp2);
                kramtrx.B3(ii+j) = 4.0 * Two_h * SSP.rho_int(j) * cs2 * (cp2 - cs2) / cp2;
                kramtrx.B4(ii+j) = Two_h * (cp2 - 2.0 * cs2) / cp2;
                kramtrx.rho(ii+j) = Two_h * std::real(omega2) * SSP.rho_int(j);
            }
        }
    }
    
    // 处理底部半空间属性
    if (params.HSBot.BC == BC_Mode::MODE_A_Half_space) {
        if (std::real(params.HSBot.cs) > 0.0) { // 弹性底部半空间
            ElasticFlag = true;
            cMin = std::min(cMin, std::real(params.HSBot.cs));
            cHigh = std::min(cHigh, std::real(params.HSBot.cs));
        } else { // 声学底部半空间
            cMin = std::min(cMin, std::real(params.HSBot.cp));
        }
    }

    // 处理顶部半空间属性
    if (params.HSTop.BC == BC_Mode::MODE_A_Half_space) {
        if (std::real(params.HSTop.cs) > 0.0) { // 弹性底部半空间
            ElasticFlag = true;
            cMin = std::min(cMin, std::real(params.HSTop.cs));
            cHigh = std::min(cHigh, std::real(params.HSTop.cs));
        } else { // 声学顶部半空间
            cMin = std::min(cMin, std::real(params.HSTop.cp));
        }
    }
    
    // 如果存在弹性介质，则减小cMin以考虑Scholte波
    if (ElasticFlag) {
        cMin = 0.85 * cMin;
    }
    Clow = std::max(Clow, cMin);
}