#ifndef KKCIFACE_H
#define KKCIFACE_H

#include "kkc_params.h"
#include <memory>

class kkc_interface_PIMPL;

// 只作为一个参数表，不涉及具体的计算过程
class kkc_interface
{
public:
    kkc_interface();     // 构造函数
    ~kkc_interface();    // 析构函数
    void init();         // 初始化
    void setup();        // 配置
    void input_setup();  // 设置输入参数
    void intm_setup();   // 设置中间矩阵参数
    void output_setup(); // 设置输出参数
    void run();          // 运行
    void clearResults(); // 清除结果
    void runField();     // 运行声场
    void runSolveV();    // 运行特征值求解器
    void free();         // 释放内存

    // 参数设置
    void set_Title(std::string &title);
    void set_Freq(double freq);
    void set_freqvec(VectorXd freqvec);           // 设置频率向量
    void set_RProf(const VectorXd &RProf);        // 设置距离剖面
    void set_RProf(const double &start, const double &end, const int &NProf); // 设置距离剖面（插值）
    void set_SSP(SSP_1D *sspInput, size_t NProf); // 设置SSP
    void set_AttenUnit(Atten_Mode mode);          // 设置衰减单位

    void set_Sz(const VectorXd &Sz);
    void set_Sz(const double &start, const double &end, const int &NSz); // 设置声源深度（插值）

    void set_Rr(const VectorXd &Rr); // 设置水平接收

    void set_Rr(const double &start, const double &end, const int &NRr);                                      // 设置水平接收（插值）
    void set_Rz(const VectorXd &Rz);                                                                          // 设置垂直接收
    void set_Rz(const double &start, const double &end, const int &NRz);                                      // 设置垂直接收（插值）
    void set_surface_Type(BC_Mode bc, size_t iprof);                                                         // 设置边界条件类型
    void set_bottom_Type(BC_Mode bc, size_t iprof);                                                          // 设置底部边界条件类型
    void set_BottomLine(double zTemp, double alphaR, double alphaI, double betaR, double betaI, double rho, size_t iprof);  // 设置底部半空间
    void set_SurfaceLine(double zTemp, double alphaR, double alphaI, double betaR, double betaI, double rho, size_t iprof); // 设置表面半空间
    void set_Clow(double Clow); // 设置最低相速度
    void set_Chigh(double Chigh); // 设置最高相速度
    void set_GridType(Grid_Mode type);                                                                        // 设置网格类型

    void set_SourceType(Source_Mode type);                          // 设置源类型
    void set_RunMode(Run_Mode mode);                                // 运行模式
    void set_Velocity_enable(bool is_Velocity);                     // 设置是否计算振速
    void set_ReflCoef_Top(std::vector<ReflectionCoef> ReflCoef);    // 设置顶部反射系数
    void set_ReflCoef_Bottom(std::vector<ReflectionCoef> ReflCoef); // 设置底部反射系数
    void set_SBP(const VectorXd &pat, const VectorXd &theta);       // 设置指向性
    VectorXd get_BottomLine(size_t iprof);                          // 获取底部半空间
    VectorXd get_SurfaceLine(size_t iprof);                         // 获取表面半空间
    double get_freq();                                              // 获取频率
    std::vector<ReflectionCoef> get_ReflCoef_Top();                 // 获取顶部反射系数
    std::vector<ReflectionCoef> get_ReflCoef_Bottom();              // 获取底部反射系数
    std::pair<VectorXd, VectorXd> get_SBP();                        // 获取指向性
    SSP_1D get_SSP(size_t iprof);                                   // 获取1D SSP
    VectorXd get_Sz();                                              // 获取声源深度
    VectorXd get_Rr();                                              // 获取水平接收
    VectorXd get_Rz();                                              // 获取垂直接收
    std::complex<float> *get_u(int srcIndex);                       // 获取某个声源复声压指针
    std::complex<float> *get_v(int srcIndex);                       // 获取某个声源复垂直振速指针
    std::complex<float> *get_h(int srcIndex);                       // 获取某个声源水平振速指针
    std::complex<float> *get_u_AllSources();                        // 获取全部声源的复声压
    std::complex<float> *get_v_AllSources();                        // 获取全部声源的垂直振速
    std::complex<float> *get_h_AllSources();                        // 获取全部声源的水平振速
    void export_result(std::string filename);                       // 导出结果到文件
    void export_shd(std::string filename, int dataType);
    parameters &getParams();                         // 获取参数的引用
    parameters getParams_Copy() const;               // 获取参数的副本
    const parameters &getParams_const() const;       // 获取参数的副本
    kkc_output &getOutput() const;                   // 获取输出的引用
    kkc_output getOutput_Copy() const;               // 获取输出的副本
    const kkc_output &getOutput_const() const;       // 获取输出的副本
    bool from_json(const std::string &jsonPath);     // 从json读取参数
    bool to_json(const std::string &jsonPath) const; // 将参数写入json
    std::string to_json_string() const;              // 将参数写入json字符串

private:
    kkc_output *output; // 输出
    TridMtx *intm_TridMtx; // 中间矩阵
    parameters *params; // 输入
    // std::atomic<int> sharedJobID;//共享任务ID

    int NumThreads;    // 线程数
    size_t totalTasks; // 任务数
    size_t field_size; // 场大小

    bool is_setup = false; // 是否设置了参数
    std::unique_ptr<kkc_interface_PIMPL> impl;
};

#endif // KKCIFACE_H
