#ifndef KKCIFACE_H
#define KKCIFACE_H

#include "OpenOceanKrakenParams.h"
#include <cstddef>
#include <memory>
#include <stdexcept>

namespace OpenOceanKraken
{
    struct FieldSnapshot
    {
        std::string title;
        double frequency = 0.0;
        Grid_Mode grid_type = Grid_Mode::MODE_R_Rectangular;
        std::size_t source_count = 0;
        std::size_t range_count = 0;
        std::size_t depth_count = 0;
        std::vector<double> source_depths;
        std::vector<double> receiver_ranges;
        std::vector<double> receiver_depths;
        std::vector<std::complex<float>> values;
    };

    struct ModeProfileSnapshot
    {
        double profile_range = 0.0;
        std::vector<double> depth;
        std::vector<std::complex<double>> wavenumbers;
        std::vector<double> group_velocity;
        Eigen::MatrixXcd mode_shapes;
    };

    class OpenOceanKraken_PIMPL;
    // 只作为一个参数表，不涉及具体的计算过程
    class KernelInterface
    {
    public:
        void set_RProf(const Eigen::VectorXd &ranges);
        void set_RProf(double start, double end, int count);
        void set_MLimit(int limit);
        void set_CoherenceType(CoherenceType type);
        void set_ModeType(ModeType type);
        void set_ModeSampling(const Eigen::VectorXd &sourceDepths,
                              const Eigen::VectorXd &receiverDepths);
        FieldSnapshot getPressureCopy() const;
        FieldSnapshot getVerticalVelocityCopy() const;
        FieldSnapshot getHorizontalVelocityCopy() const;
        std::vector<ModeProfileSnapshot> getModesCopy() const;

        KernelInterface();  // 构造函数
        KernelInterface(ThreadPool &pool); // 计算前需要传入一个线程池
        ~KernelInterface(); // 析构函数
        
        void setNumThreads(int num_threads); // 设置线程数
        void setThreadPool(ThreadPool &pool); // 设置线程池
        int getNumThreads() const; // 获取线程数
        int getHardwareThreads() const; //获取硬件线程数

        void run();                   // 运行
        void clearResults();          // 清除结果
        void runField();              // 运行声场
        void runEigen();             // 运行特征值求解器
        void free();                  // 释放内存

        // 创建参数结构
        // ssp::SSPLayer create_SSPLayer(); // 创建声速剖面层
        // ssp::SSPLayer create_SSPLayer(int npoints,
        //                               int nmesh,
        //                               double beta,
        //                               double ft,
        //                               double sigma,
        //                               const Eigen::Ref<const Eigen::VectorXd> &_z,
        //                               const Eigen::Ref<const Eigen::VectorXd> &_rho,
        //                               const Eigen::Ref<const Eigen::VectorXd> &_aR,
        //                               const Eigen::Ref<const Eigen::VectorXd> &_aI,
        //                               const Eigen::Ref<const Eigen::VectorXd> &_bR,
        //                               const Eigen::Ref<const Eigen::VectorXd> &_bI,
        //                               Media_Mode media);             // 创建声速剖面层
        // ssp::Range_Independent_Area create_Range_Independent_Area(); // 创建距离无关区域
        
        // 参数设置
        void set_Title(std::string &title);
        void set_Freq(double freq);
        void set_freqvec(Eigen::VectorXd freqvec);                              // 设置频率向量
        void set_SSP(const std::vector<ssp::Range_Independent_Area> &sspInput); // 设置SSP

        void set_AttenUnit(Atten_Mode mode); // 设置衰减单位

        void set_Sz(const Eigen::VectorXd &Sz);
        void set_Sz(const double &start, const double &end, const int &NSz); // 设置声源深度（插值）
        void set_Rr(const Eigen::VectorXd &Rr); // 设置水平接收
        void set_Rr(const double &start, const double &end, const int &NRr); // 设置水平接收（插值）
        void set_Rz(const Eigen::VectorXd &Rz);                              // 设置垂直接收
        void set_Rz(const double &start, const double &end, const int &NRz); // 设置垂直接收（插值）
        void set_Ro(const Eigen::VectorXd &Ro);                              // 设置阵列倾斜
        void set_Ro(const double &start, const double &end, const int &NRo);                              // 设置阵列倾斜
        void set_cPhase(double cLow, double cHigh);                          // 设置最低频率
        void set_GridType(Grid_Mode type);                                   // 设置网格类型
        void set_Rmax(double Rmax); // 设置最大计算距离，用于缩放error

        void set_SourceType(Source_Mode type);                                  // 设置源类型
        void set_RunMode(Run_Mode mode);                                        // 运行模式
        void set_Velocity_enable(bool is_Velocity);                             // 设置是否计算振速
        void set_ReflCoef_Top(std::vector<ReflectionCoef> ReflCoef);            // 设置顶部反射系数
        void set_ReflCoef_Bottom(std::vector<ReflectionCoef> ReflCoef);         // 设置底部反射系数
        void set_SBP(const Eigen::VectorXd &pat, const Eigen::VectorXd &theta); // 设置指向性


        std::complex<float> *get_u(int srcIndex); // 获取某个声源复声压指针
        std::complex<float> *get_v(int srcIndex); // 获取某个声源复垂直振速指针
        std::complex<float> *get_h(int srcIndex); // 获取某个声源水平振速指针
        std::complex<float> *get_u_AllSources();  // 获取全部声源的复声压
        std::complex<float> *get_v_AllSources();  // 获取全部声源的垂直振速
        std::complex<float> *get_h_AllSources();  // 获取全部声源的水平振速
        void export_result(std::string filename); // 导出振速结果到文件
        void export_mod(std::string filename);    // 导出本征值和本征函数到文件
        void export_shd(std::string filename, int dataType); // 导出声场到文件
        OOK_parameters &getParams() const;                     // 获取参数的引用
        const OOK_parameters &getParams_const() const;   // 获取参数的副本
        OOK_output &getOutput() const;                   // 获取输出的引用
        OOK_output getOutput_Copy() const;               // 获取输出的副本
        const OOK_output &getOutput_const() const;       // 获取输出的副本
#if defined(OPENOCEAN_KRAKEN_LEGACY_JSON)
        bool from_json(const std::string &jsonPath);
        bool to_json(const std::string &jsonPath) const;
        std::string to_json_string() const;
#endif
        bool from_env(const std::string &envPath);        // 从env文件读取参数 

    private:
        enum class DirtyKind
        {
            Field,
            Eigen,
            All,
            Execution
        };

        void ensureAlive() const;
        void ensureThreadPool() const;
        void ensureSetup() const;
        void releaseFieldOutput();
        void releaseEigenOutput();
        void releaseIntermediate();
        void markDirty(DirtyKind kind);
        void validateSourceIndex(int srcIndex) const;

        void init();                  // 初始化
        void setup();                 // 配置
        void input_setup();           // 设置输入参数
        void intm_setup();            // 设置中间矩阵参数
        void output_setup();          // 设置输出参数


    private:
        OOK_output *output;     // 输出
        TridMtx *intm_TridMtx;  // 中间矩阵
        OOK_parameters *params; // 输入
        ThreadPool *threadPool; // 线程池

        int NumThreads;    // 线程数
        size_t totalTasks; // 任务数
        size_t field_size; // 场大小

        bool is_setup = false; // 是否设置了参数
        
        
        std::unique_ptr<OpenOceanKraken_PIMPL> impl;
    };
}

#endif // KKCIFACE_H
