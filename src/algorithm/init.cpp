#include "init.h"



// EvaluateSSP函数声明（需要根据实际实现调整）
void EvaluateSSP(std::complex<double>* cP, std::complex<double>* cS, double* rho, 
                int Medium, int N1, double freq, const std::string& Task);

// 初始化有限差分方程
void Initialize(const parameters& params, KrakenMatrix& kramtrx) {
    bool ElasticFlag = false;
    int IAllocStat = 0;
    int Medium, NPoints, N1;

    double Two_h;
    double cP2, cS2;
    std::string Task = "TAB";
    
    // 初始化变量
    double cMin = 1e8;
    kramtrx.FirstAcoustic = 0;
    kramtrx.Loc[0] = 0; // C++使用0-based索引
    
    // 计算总网格点数
    NPoints = 0;
    kramtrx.N.resize(params.NMedia);
    for (int i = 0; i < params.NMedia; ++i) {
        NPoints += params.ssp[i].N;
        kramtrx.N(i) = params.ssp[i].N;
    }
    NPoints += params.NMedia;
    
    // 分配内存（在C++中使用vector会自动管理内存）
    kramtrx.B1.resize(NPoints);
    kramtrx.B1C.resize(NPoints);
    kramtrx.B2.resize(NPoints);
    kramtrx.B3.resize(NPoints);
    kramtrx.B4.resize(NPoints);
    kramtrx.rho.resize(NPoints);
    
    // 创建临时向量存储cP和cS
    std::vector<std::complex<double>> cP(NPoints);
    std::vector<std::complex<double>> cS(NPoints);
    
    // 处理每个介质层
    for (size_t im = 0; im < params.NMedia; ++im) { // C++使用0-based索引
        // 计算当前层的起始位置
        if (im != 0) {
            kramtrx.Loc[im] = kramtrx.Loc[im - 1] + kramtrx.N[im - 1] + 1;
        }
        
        N1 = kramtrx.N[im] + 1; // 这一层媒质的差分网格点个数
        int ii = kramtrx.Loc[im]; // C++使用0-based索引，不需要+1
        
        // 调用EvaluateSSP函数
        EvaluateSSP(&cP[ii], &cS[ii], &rho[ii], Medium + 1, N1, freq, Task); // Medium+1因为Fortran使用1-based
        
        // 加载有限差分方程的对角线
        if (std::real(cS[ii]) == 0.0) { // 声学介质情况
            params.ssp[Medium].Material = "ACOUSTIC";
            if (kramtrx.FirstAcoustic == 0) {
                kramtrx.FirstAcoustic = Medium + 1; // 保持1-based编号
            }
            kramtrx.LastAcoustic = Medium + 1; // 保持1-based编号
            
            // 计算当前层的最小声速
            double min_cP = 1e8;
            for (int j = ii; j < ii + kramtrx.N[im]; ++j) {
                min_cP = std::min(min_cP, std::real(cP[j]));
            }
            cMin = std::min(cMin, min_cP);
            
            // 计算B1和B1C
            double h_squared = kramtrx.h[im] * kramtrx.h[im];
            for (int j = ii; j < ii + kramtrx.N[im]; ++j) {
                std::complex<double> cP_squared = cP[j] * cP[j];
                std::complex<double> val = omega2 / cP_squared;
                kramtrx.B1[j] = -2.0 + h_squared * std::real(val);
                kramtrx.B1C[j] = std::imag(val);
            }
        } else { // 弹性介质情况
            if (SSP.sigma[Medium] != 0.0) {
                ERROUT("KRAKEN", "Rough elastic interfaces are not allowed");
            }
            
            SSP.Material[Medium] = "ELASTIC";
            ElasticFlag = true;
            Two_h = 2.0 * kramtrx.h[im];
            
            for (int j = ii; j < ii + kramtrx.N[im]; ++j) {
                cMin = std::min(std::real(cS[j]), cMin);
                
                cP2 = std::real(cP[j]) * std::real(cP[j]);
                cS2 = std::real(cS[j]) * std::real(cS[j]);
                
                kramtrx.B1(j) = Two_h / (kramtrx.rho[j] * cS2);
                kramtrx.B2(j) = Two_h / (kramtrx.rho[j] * cP2);
                kramtrx.B3(j) = 4.0 * Two_h * kramtrx.rho[j] * cS2 * (cP2 - cS2) / cP2;
                kramtrx.B4(j) = Two_h * (cP2 - 2.0 * cS2) / cP2;
                kramtrx.rho(j) = Two_h * std::real(omega2) * kramtrx.rho[j];
            }
        }
    }
    
    // 处理底部半空间属性
    if (HSBot.BC == "A") {
        if (std::real(HSBot.cS) > 0.0) { // 弹性底部半空间
            ElasticFlag = true;
            cMin = std::min(cMin, std::real(HSBot.cS));
            cHigh = std::min(cHigh, std::real(HSBot.cS));
        } else { // 声学底部半空间
            cMin = std::min(cMin, std::real(HSBot.cP));
            // cHigh = std::min(cHigh, std::real(HSBot.cP)); // 注释掉，与原代码一致
        }
    }
    
    // 处理顶部半空间属性
    if (HSTop.BC == "A") {
        if (std::real(HSTop.cS) > 0.0) { // 弹性顶部半空间
            ElasticFlag = true;
            cMin = std::min(cMin, std::real(HSTop.cS));
            cHigh = std::min(cHigh, std::real(HSTop.cS));
        } else { // 声学顶部半空间
            cMin = std::min(cMin, std::real(HSTop.cP));
            // cHigh = std::min(cHigh, std::real(HSTop.cP)); // 注释掉，与原代码一致
        }
    }
    
    // 如果存在弹性介质，则减小cMin以考虑Scholte波
    if (ElasticFlag) {
        cMin = 0.85 * cMin;
    }
    cLow = std::max(cLow, cMin);
}