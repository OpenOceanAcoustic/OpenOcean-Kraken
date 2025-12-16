#include "init.h"

// 初始化有限差分方程
void Initialize(int& iset, size_t iprof, parameters& params, KrakenMatrix& kramtrx, int ntimes) {
    bool ElasticFlag = false;
    int NPoints = 0;

    double Two_h;
    double cp2, cs2;
    double omega2 = SQ(2 * pi * params.freqinfo->freq);
    
    // 初始化变量
    double& Clow = params.Clow;
    double cMin = 1e8;
    double& cHigh = params.Chigh;
    params.SSP[iprof].FirstAcoustic = -1;
    params.mesh.Loc.resize(params.SSP[iprof].NMedia+1);
    params.mesh.Loc[0] = 0; // C++使用0-based索引
    
    // 计算总网格点数
    params.mesh.N.resize(params.SSP[iprof].NMedia);
    params.mesh.h.resize(params.SSP[iprof].NMedia);

    for (int i = 0; i < params.SSP[iprof].NMedia; ++i) {
        params.mesh.N(i) = params.SSP[iprof].NMesh(i) * ntimes;
        params.mesh.h(i) = params.SSP[iprof].depth(i) / params.mesh.N(i);
        NPoints += params.mesh.N(i);
        if (i == 0)
        {
            params.mesh.hV(iset) = params.mesh.h(i);
        }
    }
    NPoints += params.SSP[iprof].NMedia;
    
    // 分配内存（在C++中使用vector会自动管理内存）
    kramtrx.B1.resize(NPoints);
    kramtrx.B1C.resize(NPoints);
    kramtrx.B2.resize(NPoints);
    kramtrx.B3.resize(NPoints);
    kramtrx.B4.resize(NPoints);
    kramtrx.rho.resize(NPoints);

    SSPStructure SSP = params.SSP[iprof];
    SSP.cp_int.resize(NPoints);
    SSP.cs_int.resize(NPoints);
    SSP.rho_int.resize(NPoints);
    SSP.Material.resize(params.SSP[iprof].NMedia);
    
    // 处理每个介质层
    for (int im = 0; im < params.SSP[iprof].NMedia; ++im) { // C++使用0-based索引
        // 计算当前层的起始位置
        if (im != 0) {
            params.mesh.Loc(im) = params.mesh.Loc(im - 1) + params.mesh.N(im - 1) + 1;
        }
        
        int ii = params.mesh.Loc(im); // C++使用0-based索引，不需要+1
        Two_h = 2.0 * params.mesh.h[im];
        
        // 调用EvaluateSSP函数
        
        SSP.NMesh(im) = params.mesh.N(im);


        EvaluateSSP(SSP, params.SSP[iprof].SSPType, im);

        // 打印rho_int
        // std::cout << "RHO INT" << SSP.rho_int.transpose() << std::endl;

        
        // 加载有限差分方程的对角线
        if (std::real(SSP.cs[0]) == 0.0) { // 声学介质情况
            SSP.Material[im] = Media_Mode::MODE_A_Acoustic;
            if (params.SSP[iprof].FirstAcoustic == -1) {
                params.SSP[iprof].FirstAcoustic = im;
            }
            params.SSP[iprof].LastAcoustic = im;
            
            // 计算当前层的最小声速
            double min_cp = 1e8;
            for (int j = 0; j < SSP.NMesh(im) + 1; ++j) {
                min_cp = std::min(min_cp, std::real(SSP.cp_int(ii+j)));
            }
            cMin = std::min(cMin, min_cp);
            
            // 计算B1和B1C
            double h2 = SQ(params.mesh.h(im));
            for (int j = 0; j < SSP.NMesh(im) + 1; ++j) {
                complex<double> val = omega2 / SQ(SSP.cp_int(ii+j));
                kramtrx.B1(ii+j) = -2.0 + h2 * std::real(val);
                kramtrx.B1C(ii+j) = std::imag(val);
                kramtrx.rho(ii+j) = SSP.rho_int(ii+j);
            }

        } else { // 弹性介质情况
            SSP.Material[im] = Media_Mode::MODE_E_Elastic;
            ElasticFlag = true;
            Two_h = 2.0 * params.mesh.h[im];
            
            for (int j = 0; j < SSP.NMesh(im) + 1; ++j) {
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

    // // 打印kramtrx内部参数
    // std::cout << "B1 " << kramtrx.B1.size() << kramtrx.B1.transpose() << std::endl;

    // std::cout << "B1C " << kramtrx.B1C.size() << kramtrx.B1C.transpose() << std::endl;

    // std::cout << "rho " << kramtrx.rho.size() << kramtrx.rho.transpose() << std::endl;
    

    // 处理底部半空间属性
    if (params.HSBot[iprof].BC == BC_Mode::MODE_A_Half_space) {
        if (std::real(params.HSBot[iprof].cs) > 0.0) { // 弹性底部半空间
            ElasticFlag = true;
            cMin = std::min(cMin, std::real(params.HSBot[iprof].cs));
            cHigh = std::min(cHigh, std::real(params.HSBot[iprof].cs));
        } else { // 声学底部半空间
            cMin = std::min(cMin, std::real(params.HSBot[iprof].cp));
        }
    }

    // 处理顶部半空间属性
    if (params.HSTop[iprof].BC == BC_Mode::MODE_A_Half_space) {
        if (std::real(params.HSTop[iprof].cs) > 0.0) { // 弹性顶部半空间
            ElasticFlag = true;
            cMin = std::min(cMin, std::real(params.HSTop[iprof].cs));   
            cHigh = std::min(cHigh, std::real(params.HSTop[iprof].cs));
        } else { // 声学顶部半空间
            cMin = std::min(cMin, std::real(params.HSTop[iprof].cp));
        }
    }
    
    // 如果存在弹性介质，则减小cMin以考虑Scholte波
    if (ElasticFlag) {
        cMin = 0.85 * cMin;
    }
    Clow = std::max(Clow, cMin);
}