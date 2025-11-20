//===================================================================
// 简正波耦合模块头文件 EvaluateCMMod.h
// 
// 功能: 实现海洋声学中的简正波耦合理论，用于计算变深海洋环境中的声场
// 应用场景: 当海洋环境参数(如深度、声速剖面等)随水平距离变化时，不同位置的
//           简正波模态之间会发生能量耦合，该模块负责计算这种耦合效应下的声场
// 
// 版本: 1.0
// 作者: OpenOcean-Kraken开发团队
//===================================================================
#ifndef EVALUATECMMOD_H
#define EVALUATECMMOD_H

#include <complex>        // 用于复数运算
#include <string>         // 用于文件路径等字符串操作
#include <vector>         // 用于数据存储和处理
#include "../../include/kkc_params.h"  // 项目全局参数定义

/**
 * @class EvaluateCMMod
 * @brief 简正波耦合模块类，实现简正波耦合声场计算
 * 
 * 该类实现了海洋声学中的简正波耦合理论，用于处理水平不均匀海洋环境中的声场计算。
 * 当声线传播路径上的海洋环境参数发生变化时，不同位置的简正波模态会发生能量交换，
 * 该模块通过计算模态间的耦合系数，模拟这种复杂的声场传播现象。
 */
class EvaluateCMMod {
public:
    /**
     * @brief 构造函数
     * 
     * 初始化简正波耦合模块，设置内部状态变量的初始值
     */
    EvaluateCMMod();
    
    /**
     * @brief 析构函数
     * 
     * 清理模块资源，释放可能的动态分配内存
     */
    ~EvaluateCMMod();
    
    /**
     * @brief 核心函数：使用简正波耦合方法计算压力场
     * 
     * 该函数是模块的主要入口点，通过简正波耦合理论计算水平不均匀海洋环境中的声场分布
     * 
     * @param fileRoot 数据文件根路径，用于读取模态信息
     * @param rProf 不同水平位置处的距离数组
     * @param nProf 水平位置的数量
     * @param phiS 源位置的模态函数
     * @param phi 存储各位置模态函数的二维数组
     * @param rd 接收深度数组
     * @param nrd 接收深度数量
     * @param r 计算声场的距离数组
     * @param nr 计算声场的距离点数量
     * @param k 波数数组
     * @param m 模态数量
     * @param option 计算选项标志
     * @param p 输出参数：计算得到的压力场
     */
    void evaluateCM(const std::string& fileRoot, std::vector<double>& rProf, const int& nProf, 
                    const std::vector<std::complex<double>>& phiS, std::vector<std::vector<std::complex<double>>>& phi, 
                    const std::vector<double>& rd, const int& nrd, const std::vector<double>& r, const int& nr, 
                    std::vector<std::complex<double>>& k, int& m, const std::string& option, 
                    std::vector<std::vector<std::complex<double>>>& p);

private:
    // 私有常量
    const int MAX_M = 20000;     ///< 最大模态数量
    const int MAX_N = 100001;    ///< 最大深度采样点数
    const int MAX_MED = 501;     ///< 最大介质层数
    const int MAX_NFREQ = 1000;  ///< 最大频率数量
    const double PI = 3.1415926; ///< 圆周率近似值
    
    // 私有成员变量
    bool firstCall;     ///< 首次调用标志，用于初始化操作
    int iRecProfileR;   ///< 文件记录指针，用于模态文件读取
    
    // 私有辅助函数
    /**
     * @brief 处理新剖面对应的模态信息
     * 
     * 当处理到新的水平位置时，读取该位置的模态信息，并计算模态振幅系数
     * 
     * @param fileRoot 数据文件根路径
     * @param ifreq 当前频率索引
     * @param k 波数数组，输出参数
     * @param phiR 模态函数矩阵，输出参数
     * @param mr 模态数量，输出参数
     * @param rd 接收深度数组
     * @param nrd 接收深度数量
     * @param a 模态振幅系数，输出参数
     */
    void newProfile(const std::string& fileRoot, const int& ifreq, std::vector<std::complex<double>>& k, 
                    std::vector<std::vector<std::complex<double>>>& phiR, int& mr, 
                    const std::vector<double>& rd, const int& nrd, std::vector<std::complex<double>>& a);
    
    /**
     * @brief 计算界面左侧的压力场
     * 
     * 在两个水平界面之间的左侧位置，计算该位置处的声场分布
     * 
     * @param fileRoot 数据文件根路径
     * @param ifreq 当前频率索引
     * @param iRecProfile 记录指针
     * @param a 模态振幅系数数组
     * @param k 波数数组
     * @param z 深度数组
     * @param m 模态数量
     * @param p 压力场，输出参数
     * @param nr 接收点数量
     * @param nTot 总深度点数
     * @param bcTop 顶部边界条件
     * @param rhoT 顶部密度
     * @param kTop2 顶部波数平方
     * @param depthT 顶部深度
     * @param bcBot 底部边界条件
     * @param rhoB 底部密度
     * @param kBot2 底部波数平方
     * @param depthB 底部深度
     * @param gamTL 顶部水平波数
     * @param gamBL 底部水平波数
     * @param depthTL 顶部有效深度
     * @param depthBL 底部有效深度
     * @param phiTL 顶部模态函数
     * @param phiBL 底部模态函数
     * @param ml 局部模态数量
     */
    void pLeft(const std::string& fileRoot, const int& ifreq, int& iRecProfile, 
               std::vector<std::complex<double>>& a, std::vector<std::complex<double>>& k, 
               std::vector<double>& z, int& m, std::vector<std::complex<double>>& p, 
               int& nr, int& nTot, char& bcTop, double& rhoT, std::complex<double>& kTop2, 
               double& depthT, char& bcBot, double& rhoB, std::complex<double>& kBot2, 
               double& depthB, std::vector<std::complex<double>>& gamTL, 
               std::vector<std::complex<double>>& gamBL, double& depthTL, double& depthBL, 
               std::vector<std::complex<double>>& phiTL, std::vector<std::complex<double>>& phiBL, 
               int& ml);
    
    /**
     * @brief 计算半空间尾部贡献
     * 
     * 计算半空间区域对总声场的尾部贡献，用于处理非理想边界条件
     * 
     * @param depth 计算点深度
     * @param phiL 左侧模态函数数组
     * @param gamL 左侧衰减常数数组
     * @param depthL 左侧深度
     * @param ml 局部模态数量
     * @param phiR 右侧模态函数值
     * @param gamR 右侧衰减常数值
     * @param depthR 右侧深度
     * @return 尾部贡献的复数振幅值
     */
    std::complex<double> calculateTail(const double& depth, const std::vector<std::complex<double>>& phiL, 
                                       const std::vector<std::complex<double>>& gamL, const double& depthL, 
                                       const int& ml, const std::complex<double>& phiR, 
                                       const std::complex<double>& gamR, const double& depthR);
    
    /**
     * @brief 求解Pekeris波导的特征根
     * 
     * 计算Pekeris波导中给定gamma^2值对应的特征根，用于特定波导模型的模态计算
     * 
     * @param gamma2 gamma的平方值
     * @return 计算得到的特征根
     */
    std::complex<double> pekerisRoot(const std::complex<double>& gamma2);
    
    /**
     * @brief 从文件中读取模态信息
     * 
     * 读取指定位置和频率的模态信息，包括波数和模态函数
     * 
     * @param fileRoot 数据文件根路径
     * @param iProf 位置索引
     * @param ifreq 频率索引
     * @param maxM 最大模态数量
     * @param rd 接收深度数组
     * @param nrd 接收深度数量
     * @param option 计算选项标志
     * @param k 波数数组，输出参数
     * @param phi 模态函数矩阵，输出参数
     * @param m1 模态数量，输出参数
     * @param freqVec 频率数组
     * @param nfreq 频率数量
     * @param title 标题信息
     */
    void getModes(const std::string& fileRoot, const int& iProf, const int& ifreq, const int& maxM, 
                 const std::vector<double>& rd, const int& nrd, const std::string& option, 
                 std::vector<std::complex<double>>& k, std::vector<std::vector<std::complex<double>>>& phi, 
                 int& m1, std::vector<double>& freqVec, int& nfreq, std::string& title);
    
    /**
     * @brief 读取模态文件头信息
     * 
     * 读取模态文件的头部信息，包括频率、介质参数等
     * 
     * @param fileRoot 数据文件根路径
     * @param iProf 位置索引
     * @param iRecProfile 记录指针，输出参数
     * @param lRecl 记录长度，输出参数
     * @param title 标题信息，输出参数
     * @param freqVec 频率数组，输出参数
     * @param nfreq 频率数量，输出参数
     * @param nMedia 介质数量，输出参数
     * @param nl 层数，输出参数
     * @param nMat 材料数量，输出参数
     * @param n 各层采样点数，输出参数
     * @param material 材料名称数组，输出参数
     * @param depth 层深度数组，输出参数
     * @param rho 密度数组，输出参数
     * @param z 深度坐标数组，输出参数
     */
    void readModeHeader(const std::string& fileRoot, const int& iProf, int& iRecProfile, 
                       int& lRecl, std::string& title, std::vector<double>& freqVec, 
                       int& nfreq, int& nMedia, int& nl, int& nMat, 
                       std::vector<int>& n, std::vector<std::string>& material, 
                       std::vector<double>& depth, std::vector<double>& rho, 
                       std::vector<double>& z);
    
    /**
     * @brief 读取波数信息
     * 
     * 从模态文件中读取指定频率的波数信息
     * 
     * @param iRecProfile 记录指针
     * @param ifreq 频率索引
     * @param k 波数数组，输出参数
     * @param ml 模态数量，输出参数
     * @param maxM 最大模态数量
     * @param lRecl 记录长度
     */
    void readWavenumbers(const int& iRecProfile, const int& ifreq, 
                        std::vector<std::complex<double>>& k, int& ml, const int& maxM, 
                        const int& lRecl);
    
    /**
     * @brief 计算积分权重
     * 
     * 计算在接收深度上的积分权重，用于数值积分计算
     * 
     * @param z 深度数组
     * @param nTot 总深度点数
     * @param rd 接收深度数组
     * @param nrd 接收深度数量
     * @param w 权重数组，输出参数
     * @param ird 接收深度索引数组，输出参数
     */
    void weight(const std::vector<double>& z, const int& nTot, const std::vector<double>& rd, 
               const int& nrd, std::vector<double>& w, std::vector<int>& ird);
};

#endif // EVALUATECMMOD_H