// 简正波耦合模块头文件
// 用于海洋声学计算中的简正波耦合理论实现
#ifndef EVALUATECMMOD_H
#define EVALUATECMMOD_H

#include <complex>
#include <string>
#include <vector>
#include "../../include/kkc_params.h"

class EvaluateCMMod {
public:
    // 构造函数和析构函数
    EvaluateCMMod();
    ~EvaluateCMMod();
    
    // 主要函数：计算压力场
    void evaluateCM(const std::string& fileRoot, std::vector<double>& rProf, const int& nProf, 
                    const std::vector<std::complex<double>>& phiS, std::vector<std::vector<std::complex<double>>>& phi, 
                    const std::vector<double>& rd, const int& nrd, const std::vector<double>& r, const int& nr, 
                    std::vector<std::complex<double>>& k, int& m, const std::string& option, 
                    std::vector<std::vector<std::complex<double>>>& p);

private:
    // 私有常量
    const int MAX_M = 20000;
    const int MAX_N = 100001;
    const int MAX_MED = 501;
    const int MAX_NFREQ = 1000;
    const double PI = 3.1415926;
    
    // 私有成员变量
    bool firstCall;  // 首次调用标志
    int iRecProfileR; // 记录指针
    
    // 私有辅助函数
    void newProfile(const std::string& fileRoot, const int& ifreq, std::vector<std::complex<double>>& k, 
                    std::vector<std::vector<std::complex<double>>>& phiR, int& mr, 
                    const std::vector<double>& rd, const int& nrd, std::vector<std::complex<double>>& a);
    
    void pLeft(const std::string& fileRoot, const int& ifreq, int& iRecProfile, 
               std::vector<std::complex<double>>& a, std::vector<std::complex<double>>& k, 
               std::vector<double>& z, int& m, std::vector<std::complex<double>>& p, 
               int& nr, int& nTot, char& bcTop, double& rhoT, std::complex<double>& kTop2, 
               double& depthT, char& bcBot, double& rhoB, std::complex<double>& kBot2, 
               double& depthB, std::vector<std::complex<double>>& gamTL, 
               std::vector<std::complex<double>>& gamBL, double& depthTL, double& depthBL, 
               std::vector<std::complex<double>>& phiTL, std::vector<std::complex<double>>& phiBL, 
               int& ml);
    
    std::complex<double> calculateTail(const double& depth, const std::vector<std::complex<double>>& phiL, 
                                       const std::vector<std::complex<double>>& gamL, const double& depthL, 
                                       const int& ml, const std::complex<double>& phiR, 
                                       const std::complex<double>& gamR, const double& depthR);
    
    // 辅助工具函数
    std::complex<double> pekerisRoot(const std::complex<double>& gamma2);
    
    // 用于从文件读取模态信息的辅助函数
    void getModes(const std::string& fileRoot, const int& iProf, const int& ifreq, const int& maxM, 
                 const std::vector<double>& rd, const int& nrd, const std::string& option, 
                 std::vector<std::complex<double>>& k, std::vector<std::vector<std::complex<double>>>& phi, 
                 int& m1, std::vector<double>& freqVec, int& nfreq, std::string& title);
    
    void readModeHeader(const std::string& fileRoot, const int& iProf, int& iRecProfile, 
                       int& lRecl, std::string& title, std::vector<double>& freqVec, 
                       int& nfreq, int& nMedia, int& nl, int& nMat, 
                       std::vector<int>& n, std::vector<std::string>& material, 
                       std::vector<double>& depth, std::vector<double>& rho, 
                       std::vector<double>& z);
    
    void readWavenumbers(const int& iRecProfile, const int& ifreq, 
                        std::vector<std::complex<double>>& k, int& ml, const int& maxM, 
                        const int& lRecl);
    
    void weight(const std::vector<double>& z, const int& nTot, const std::vector<double>& rd, 
               const int& nrd, std::vector<double>& w, std::vector<int>& ird);
};

#endif // EVALUATECMMOD_H