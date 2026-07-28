#include "OpenOceanKrakenInterface.h"
#include "util.h"
#include "input_SSP.hpp"
#include "input_Sz_Rz_RR.hpp"
#include "input_reflcoef.hpp"
#include "input_sbp.hpp"
#include "input_Freq.hpp"
#include "output_eigen.hpp"
#include "output_field.hpp"
#include "run.h"
#include "EvaluateAD.h"
#include "EvaluateCM.h"
#include "json_in_out.hpp"
#include "env_in_out.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstring>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace OpenOceanKraken
{
    void print_header()
    {
        std::cout << " - OpenOcean-Kraken: A parallel underwater acoustic nomal modes simulator with pybind11 Python and C++ interfaces\n"
                     "\n"
                     "Copyright (C) 2024-2025 OpenOceanAcoustic\n"
                     "Authors: Qian Peng,Liang Yi, Dai Nuoge, Luo Shixiong, Liu Ruihang, Liu Jiatong, Wang Zhengwei\n"
                     "Mail: yingxinliang1@gmail.com\n"
                     "Based on Kraken, which is Copyright (C) 1983-2024 Michael B. Porter\n"
                     "GPL3 licensed, no warranty, see LICENSE or https://www.gnu.org/licenses/\n"
                     "\n";
    }
    // 定义实现类
    class OpenOceanKraken_PIMPL
    {
    public:
        input_Freq INPUT_FREQ;
        input_SSP INPUT_SSP;
        input_Sz_Rz_RR INPUT_SZ_RZ_RR;
        input_reflcoef INPUT_REFLCOEF;
        input_sbp INPUT_SBP;

        output_Eigen OUTPUT_EIGEN;
        output_Field OUTPUT_FIELD;
        OpenOceanKraken_PIMPL()
        {
        }
        ~OpenOceanKraken_PIMPL()
        {
        }
    };

    Interface::Interface()
        : params(new OOK_parameters()),
          output(new OOK_output()),
          intm_TridMtx(nullptr),  // 初始化为 nullptr
          threadPool(nullptr),
          NumThreads(1),
          totalTasks(0),
          field_size(0),
          impl(std::make_unique<OpenOceanKraken_PIMPL>())
    {
        this->init(); // 初始化
    }
    Interface::Interface(ThreadPool &pool)
        : params(new OOK_parameters()),
          output(new OOK_output()),
          intm_TridMtx(nullptr),  // 初始化为 nullptr
          NumThreads(1),
          totalTasks(0),
          field_size(0),
          impl(std::make_unique<OpenOceanKraken_PIMPL>()),
          threadPool(&pool) // 初始化传入的线程池引用
    {
        this->init(); // 初始化
    }

    Interface::~Interface()
    {
        // 释放params的各个成员变量的内存
        this->free(); // 释放内存
    }

    void Interface::init()
    {
        if (!this->params)
        {
            std::cerr << "Error: Failed to allocate memory for parameters" << std::endl;
            exit(1);
        }
        this->field_size = 0; // 场大小初始化为0

        this->is_setup = false;           // 初始化是否设置参数标志位
        auto &params = this->getParams(); // 获取参数
        auto &output = this->getOutput(); // 获取输出

        // 参数初始化
        impl->INPUT_FREQ.Init(params);
        impl->INPUT_SSP.Init(params);
        impl->INPUT_SZ_RZ_RR.Init(params);
        impl->INPUT_REFLCOEF.Init(params);
        impl->INPUT_SBP.Init(params);

        // default params
        impl->INPUT_FREQ.Default(params);
        impl->INPUT_SSP.Default(params);
        impl->INPUT_SZ_RZ_RR.Default(params);
        impl->INPUT_REFLCOEF.Default(params);
        impl->INPUT_SBP.Default(params);

        // output初始化
        impl->OUTPUT_FIELD.Init(output);
        impl->OUTPUT_EIGEN.Init(output);
        
    }

    void Interface::input_setup()
    {
        ensureAlive();
        auto &params = this->getParams(); // 获取参数
        auto &output = this->getOutput(); // 获取输出
        // 预处理
        impl->INPUT_SSP.Preprocess(params);
        impl->INPUT_SZ_RZ_RR.Preprocess(params);
        impl->INPUT_REFLCOEF.Preprocess(params);
        impl->INPUT_SBP.Preprocess(params);
        impl->INPUT_FREQ.Preprocess(params);

        // 初始化声压共享内存区域
        this->field_size = (size_t)params.Pos.NSz * (size_t)params.Pos.NRz_per_range * (size_t)params.Pos.NRr;
        std::cout << "OpenOcean-Kraken: input field size: " << this->field_size << std::endl;
    }
    void Interface::intm_setup()
    {
        ensureAlive();
        auto &params = this->getParams(); // 获取参数
        size_t workspace_count = std::max<size_t>(1, params.SSP.size());
        if (this->NumThreads == 0)
        {
            throw std::runtime_error("NumThreads must be greater than zero.");
        }
        releaseIntermediate();
        this->intm_TridMtx = new TridMtx[workspace_count];
        for (size_t i = 0; i < workspace_count; i++)
        {
            this->intm_TridMtx[i].resize(params.NMeshMax, params.NMediaMax, params.mesh.NSets);
        }
        std::cout << "OpenOcean-Kraken: intm TridMtx  " << std::endl;
    }

    // 输出配置
    void Interface::output_setup() // 设置输出参数
    {
        auto &params = this->getParams(); // 获取参数
        auto &output = this->getOutput(); // 获取输出
        impl->OUTPUT_EIGEN.Preprocess(params, output);
        impl->OUTPUT_FIELD.Preprocess(params, output);
    }

    void Interface::ensureAlive() const
    {
        if (!this->params || !this->output)
        {
            throw std::runtime_error("OpenOceanKraken::Interface has already been freed.");
        }
    }

    void Interface::ensureThreadPool() const
    {
        ensureAlive();
        if (!this->threadPool)
        {
            throw std::runtime_error("ThreadPool is not set. Construct Interface with ThreadPool or call setThreadPool() before run().");
        }
        if (this->NumThreads <= 0)
        {
            throw std::runtime_error("NumThreads must be greater than zero.");
        }
    }

    void Interface::ensureSetup() const
    {
        ensureAlive();
        if (!this->is_setup || !this->intm_TridMtx || !this->output->eigen || !this->output->u_AllSources)
        {
            throw std::runtime_error("Interface output is not ready. Call run() before reading field output.");
        }
    }

    void Interface::releaseFieldOutput()
    {
        if (this->output)
        {
            this->impl->OUTPUT_FIELD.Finalize(*this->output);
        }
    }

    void Interface::releaseEigenOutput()
    {
        if (this->output)
        {
            this->impl->OUTPUT_EIGEN.Finalize(*this->output);
        }
    }

    void Interface::releaseIntermediate()
    {
        delete[] this->intm_TridMtx;
        this->intm_TridMtx = nullptr;
    }

    void Interface::markDirty(DirtyKind kind)
    {
        if (!this->output)
        {
            this->is_setup = false;
            return;
        }

        switch (kind)
        {
        case DirtyKind::Field:
            releaseFieldOutput();
            break;
        case DirtyKind::Eigen:
            releaseFieldOutput();
            releaseEigenOutput();
            releaseIntermediate();
            break;
        case DirtyKind::All:
            releaseFieldOutput();
            releaseEigenOutput();
            releaseIntermediate();
            break;
        case DirtyKind::Execution:
            releaseIntermediate();
            break;
        }
        this->is_setup = false;
        this->field_size = 0;
    }

    void Interface::validateSourceIndex(int srcIndex) const
    {
        ensureSetup();
        const auto &input = this->getParams_const();
        if (srcIndex < 0 || srcIndex >= input.Pos.NSz)
        {
            std::ostringstream oss;
            oss << "Source index " << srcIndex << " is out of range [0, " << input.Pos.NSz << ").";
            throw std::out_of_range(oss.str());
        }
    }

    // 配置
    void Interface::setup() // 设置参数
    {
        input_setup(); // 设置输入参数(外界输入参数，因此每次计算时候都重新计算一次环境)
        intm_setup();  // 设置中间矩阵参数

        if (!this->is_setup)
        {
            this->output_setup();  // 若未配置，则配置输出的内存
            this->is_setup = true; // 设置配置标志位
        }
    }

    void Interface::setNumThreads(int num_threads)
    {
        ensureAlive();
        if (num_threads <= 0)
        {
            throw std::invalid_argument("NumThreads must be greater than zero.");
        }
        if (this->NumThreads != num_threads)
        {
            markDirty(DirtyKind::Execution);
        }
        this->NumThreads = num_threads;
    }
    void Interface::setThreadPool(ThreadPool &pool)
    {
        ensureAlive();
        this->threadPool = &pool;
    }

    int Interface::getNumThreads() const
    {
        ensureAlive();
        return this->NumThreads;
    }

    int Interface::getHardwareThreads() const
    {
        return std::thread::hardware_concurrency();
    }

    void Interface::clearResults() // 清除结果
    {
        auto &output = this->getOutput(); // 获取输出
        auto &params = this->getParams(); // 获取参数

        this->impl->OUTPUT_FIELD.ClearResults(params, output);
        this->impl->OUTPUT_EIGEN.ClearResults(params, output);
    }
    void Interface::runField() // 运行声场
    {
        // FieldSolveWorker();
        ensureAlive();
        if (!this->is_setup)
        {
            this->setup();
        }
        auto &params = this->getParams(); // 获取参数
        auto &output = this->getOutput(); // 获取输出
        if (params.runMode == Run_Mode::MODE_M_Modes)
        {
            return;
        }
        if (params.SSP.size() > 1)
        {
            if (params.is_Velocity)
            {
                throw std::logic_error("multi-profile velocity is not implemented");
            }
            for (int isz = 0; isz < params.Pos.NSz; ++isz)
            {
                if (params.modeType == ModeType::Adiabatic)
                {
                    EvaluateAD(output.eigen, params.NProf, params, isz, output.u_AllSources);
                }
                else if (params.modeType == ModeType::Couple)
                {
                    EvaluateCM(output.eigen, params.NProf, params, isz, output.u_AllSources);
                }
                else
                {
                    throw std::logic_error("unsupported multi-profile field mode");
                }
            }
            return;
        }
        for (size_t iprof = 0; iprof < params.SSP.size(); iprof++)
        {
            FieldWorker(iprof, params, output);
        }
    }

    void Interface::runEigen() // 运行特征值
    {
        auto &params = this->getParams(); // 获取参数
        auto &output = this->getOutput(); // 获取输出
        ensureThreadPool();
        if (!this->is_setup)
        {
            this->setup();
        }
        this->clearResults();
        for (size_t iprof = 0; iprof < params.SSP.size(); iprof++)
        {
            EigenVWorker(*this->threadPool, this->NumThreads, iprof, params, this->intm_TridMtx[iprof], output);
        }
    }

    void Interface::run() // 运行
    {
        print_header();                                                                                           // 打印头信息
        this->setup();                                                                                            // 配置
        auto start = std::chrono::high_resolution_clock::now();                                                   // 计时开始
        this->runEigen();                                                                                         // 求解本征值和本征函数
        auto end = std::chrono::high_resolution_clock::now();                                                     // 计时结束
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);                       // 计算时间差
        std::cout << "Eigenvalue and Eigenfunction calculation time: " << duration.count() << " ms" << std::endl; // 打印时间差
        start = std::chrono::high_resolution_clock::now();                                                        // 计时开始
        this->runField();                                                                                         // 求解声场
        end = std::chrono::high_resolution_clock::now();                                                          // 计时结束
        duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);                            // 计算时间差
        std::cout << "Field calculation time: " << duration.count() << " ms" << std::endl;                        // 打印时间差
    }

    void Interface::free() // 释放内存
    {
        // 释放output的各个成员变量的内存
        if (!this->params && !this->output && !this->intm_TridMtx)
        {
            return;
        }
        releaseFieldOutput();
        releaseEigenOutput();
        releaseIntermediate();
        delete this->params;
        delete this->output;
        this->params = nullptr;
        this->output = nullptr;
        this->threadPool = nullptr;
        this->is_setup = false;
        this->field_size = 0;
    }

    // 参数设置
    void Interface::set_Title(std::string &title)
    {
        auto &params = this->getParams(); // 获取参数
        this->impl->INPUT_FREQ.set_title(params, title);
        markDirty(DirtyKind::All);
    }
    void Interface::set_Freq(double freq)
    {
        auto &params = this->getParams(); // 获取参数
        this->impl->INPUT_FREQ.set_freq(params, freq);
        markDirty(DirtyKind::All);
    }
    void Interface::set_freqvec(Eigen::VectorXd freqvec) // 设置频率向量
    {
        auto &params = this->getParams(); // 获取参数
        this->impl->INPUT_FREQ.set_freqvec(params, freqvec);
        markDirty(DirtyKind::All);
    }
    void Interface::set_SSP(const std::vector<ssp::Range_Independent_Area> &sspInput) // 设置SSP
    {
        auto &params = this->getParams(); // 获取参数
        this->impl->INPUT_SSP.set_SSP(params, sspInput);
        markDirty(DirtyKind::All);
    }

    void Interface::set_AttenUnit(Atten_Mode mode) // 设置衰减单位
    {
        auto &params = this->getParams(); // 获取参数
        this->impl->INPUT_SSP.set_AttenUnit(params, mode);
        markDirty(DirtyKind::All);
    }

    void Interface::set_Sz(const Eigen::VectorXd &Sz)
    {
        auto &params = this->getParams(); // 获取参数
        this->impl->INPUT_SZ_RZ_RR.set_Sz(params, Sz);
        markDirty(DirtyKind::All);
    }
    void Interface::set_Sz(const double &start, const double &end, const int &NSz) // 设置声源深度（插值）
    {
        auto &params = this->getParams(); // 获取参数
        this->impl->INPUT_SZ_RZ_RR.set_Sz(params, start, end, NSz);
        markDirty(DirtyKind::All);
    }
    void Interface::set_Rr(const Eigen::VectorXd &Rr) // 设置水平接收
    {
        auto &params = this->getParams(); // 获取参数
        this->impl->INPUT_SZ_RZ_RR.set_Rr(params, Rr);
        markDirty(DirtyKind::All);
    }
    void Interface::set_Rr(const double &start, const double &end, const int &NRr) // 设置水平接收（插值）
    {
        auto &params = this->getParams(); // 获取参数
        this->impl->INPUT_SZ_RZ_RR.set_Rr(params, start, end, NRr);
        markDirty(DirtyKind::All);
    }
    void Interface::set_Rz(const Eigen::VectorXd &Rz)
    {
        auto &params = this->getParams(); // 获取参数
        this->impl->INPUT_SZ_RZ_RR.set_Rz(params, Rz);
        markDirty(DirtyKind::All);
    }
    void Interface::set_Rz(const double &start, const double &end, const int &NRz) // 设置接收深度（插值）
    {
        auto &params = this->getParams(); // 获取参数
        this->impl->INPUT_SZ_RZ_RR.set_Rz(params, start, end, NRz);
        markDirty(DirtyKind::All);
    }
    void Interface::set_Ro(const Eigen::VectorXd &Ro) // 设置阵列倾斜
    {
        auto &params = this->getParams(); // 获取参数
        this->impl->INPUT_SZ_RZ_RR.set_Ro(params, Ro);
        markDirty(DirtyKind::All);
    }
    void Interface::set_Ro(const double &start, const double &end, const int &NRo) // 设置阵列倾斜
    {
        auto &params = this->getParams(); // 获取参数
        this->impl->INPUT_SZ_RZ_RR.set_Ro(params, start, end, NRo);
        markDirty(DirtyKind::All);
    }

    void Interface::set_cPhase(double cLow, double cHigh) // 设置最低频率
    {
        auto &params = this->getParams(); // 获取参数
        this->impl->INPUT_FREQ.set_cPhase(params, cLow, cHigh);
        markDirty(DirtyKind::All);
    }
    void Interface::set_GridType(Grid_Mode type) // 设置网格类型
    {
        auto &params = this->getParams(); // 获取参数
        this->impl->INPUT_SZ_RZ_RR.set_GridType(params, type);
        markDirty(DirtyKind::All);
    }

    void Interface::set_Rmax(double Rmax) // 设置最大计算距离，用于缩放error
    {
        auto &params = this->getParams(); // 获取参数
        this->impl->INPUT_SZ_RZ_RR.set_RMax(params, Rmax);
        markDirty(DirtyKind::All);
    }
    void Interface::set_SourceType(Source_Mode type) // 设置源类型
    {
        auto &params = this->getParams(); // 获取参数
        this->impl->INPUT_FREQ.set_SourceType(params, type);
        markDirty(DirtyKind::All);
    }
    void Interface::set_RunMode(Run_Mode mode)
    {
        auto &params = this->getParams(); // 获取参数
        this->impl->INPUT_FREQ.set_RunMode(params, mode);
        markDirty(DirtyKind::Field);
    } // 运行模式
    void Interface::set_Velocity_enable(bool is_Velocity) // 设置是否计算振速
    {
        auto &params = this->getParams(); // 获取参数
        this->impl->INPUT_FREQ.set_Velocity_enable(params, is_Velocity);
        markDirty(DirtyKind::Field);
    }
    void Interface::set_ReflCoef_Top(std::vector<ReflectionCoef> ReflCoef)
    {
        auto &params = this->getParams(); // 获取参数
        this->impl->INPUT_REFLCOEF.set_ReflCoef_Top(params, ReflCoef);
        markDirty(DirtyKind::All);
    } // 设置顶部反射系数
    void Interface::set_ReflCoef_Bottom(std::vector<ReflectionCoef> ReflCoef)
    {
        auto &params = this->getParams(); // 获取参数
        this->impl->INPUT_REFLCOEF.set_ReflCoef_Bot(params, ReflCoef);
        markDirty(DirtyKind::All);
    } // 设置底部反射系数
    void Interface::set_SBP(const Eigen::VectorXd &pat, const Eigen::VectorXd &theta) // 设置指向性
    {
        auto &params = this->getParams(); // 获取参数
        this->impl->INPUT_SBP.set_Pat(params, pat, theta);
        markDirty(DirtyKind::All);
    }

    bool Interface::to_json(const std::string &jsonPath) const // 将参数写入json
    {
        
        auto &params = this->getParams_const(); // 获取参数
        std::string jsonStr = this->to_json_string();
        std::ofstream ofs(jsonPath);
        if (!ofs.is_open())
        {
            std::cerr << "Error: Failed to open file " << jsonPath << std::endl;
            return false;
        }
        ofs << jsonStr << std::endl;
        ofs.close();

        return true;

    }
    bool Interface::from_json(const std::string &jsonPath) // 从json读取参数
    {

        std::ifstream ifs(jsonPath);
        if (!ifs.is_open())
        {
            std::cerr << "Error: Failed to open JSON file " << jsonPath << std::endl;
            return false;
        }
        try
        {
            OpenOcean_json json;
            ifs >> json;
            OOK_parameters candidate = json.get<OOK_parameters>();
            this->impl->INPUT_SSP.set_SSP(
                candidate, candidate.sspInput);
            this->getParams() = std::move(candidate);
            markDirty(DirtyKind::All);
            return true;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Error: Failed to parse JSON file " << jsonPath << ": " << e.what() << std::endl;
            return false;
        }
    }


    std::string Interface::to_json_string() const // 将参数写入json字符串
    {
        auto &params = this->getParams_const(); // 获取参数
        OpenOcean_json json;
        json = params;
        return json.dump(4);
    }

    bool Interface::from_env(const std::string &envPath)
    {
        try
        {
            OOK_parameters candidate;
            if (!read_env_file(envPath, candidate))
                return false;
            if (!read_flp_file(envPath, candidate))
                return false;

            this->impl->INPUT_SSP.set_SSP(
                candidate, candidate.sspInput);
            this->getParams() = std::move(candidate);
            markDirty(DirtyKind::All);
            return true;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Error: Failed to parse ENV file "
                      << envPath << ": " << e.what() << std::endl;
            return false;
        }
    }

    std::complex<float> *Interface::get_u(int srcIndex) // 获取某个声源复声压指针
    {
        auto &output = this->getOutput_const(); // 获取输出
        auto &input = this->getParams_const();  // 获取输入
        validateSourceIndex(srcIndex);
        if (!output.u_AllSources)
        {
            throw std::runtime_error("Pressure output is not available.");
        }
        return output.u_AllSources + srcIndex * input.Pos.NRr * input.Pos.NRz_per_range;
    }

    std::complex<float> *Interface::get_v(int srcIndex) // 获取某个声源复垂直振速
    {
        auto &output = this->getOutput_const(); // 获取输出
        auto &input = this->getParams_const();  // 获取输入
        validateSourceIndex(srcIndex);
        if (!output.v_AllSources)
        {
            throw std::runtime_error("Vertical velocity output is not available.");
        }
        return output.v_AllSources + srcIndex * input.Pos.NRr * input.Pos.NRz_per_range;
    }

    std::complex<float> *Interface::get_h(int srcIndex) // 获取某个声源水平振速
    {
        auto &output = this->getOutput_const(); // 获取输出
        auto &input = this->getParams_const();  // 获取输入
        validateSourceIndex(srcIndex);
        if (!output.h_AllSources)
        {
            throw std::runtime_error("Horizontal velocity output is not available.");
        }
        return output.h_AllSources + srcIndex * input.Pos.NRr * input.Pos.NRz_per_range;
    }

    std::complex<float> *Interface::get_u_AllSources() // 获取全部声源的复声压
    {
        auto &output = this->getOutput_const(); // 获取输出
        ensureSetup();
        if (!output.u_AllSources)
        {
            throw std::runtime_error("Pressure output is not available.");
        }
        return output.u_AllSources;
    }

    std::complex<float> *Interface::get_v_AllSources() // 获取全部声源的垂直振速
    {
        auto &output = this->getOutput_const(); // 获取输出
        ensureSetup();
        if (!output.v_AllSources)
        {
            throw std::runtime_error("Vertical velocity output is not available.");
        }
        return output.v_AllSources;
    }

    std::complex<float> *Interface::get_h_AllSources() // 获取全部声源的水平振速
    {
        auto &output = this->getOutput_const(); // 获取输出
        ensureSetup();
        if (!output.h_AllSources)
        {
            throw std::runtime_error("Horizontal velocity output is not available.");
        }
        return output.h_AllSources;
    }

    void Interface::export_result(std::string filename) // 导出结果到文件
    {
        // std::cout << "开始导出结果..." << std::endl;
        auto &input = this->getParams_const(); // 获取输入

        // TODO 导出mod
        // std::cout << "导出声压场数据..." << std::endl;
        std::string filename_P = filename + "_P";
        // this->export_Pressure(filename_P); // 导出声压场数据
        this->export_mod(filename_P);
        this->export_shd(filename_P, 1);
        if (input.is_Velocity)
        {
            std::string filename_V = filename + "_V";
            std::string filename_H = filename + "_H";
            this->export_shd(filename_V, 2); // 导出垂直振速
            this->export_shd(filename_H, 3); // 导出水平振速
        }
    }

    void Interface::export_mod(std::string filename) // 导出模型到文件
    {
        auto &input = this->getParams_const();
        auto &output = this->getOutput_const();
        if (!output.eigen)
        {
            throw std::runtime_error("Cannot export MOD file before eigen output is available.");
        }
        if (!this->intm_TridMtx)
        {
            throw std::runtime_error("Cannot export MOD file before the tridiagonal workspace is available.");
        }
        if (input.SSP.empty())
        {
            throw std::runtime_error("Cannot export MOD file without SSP profiles.");
        }

        const int nfreq = 1; // OOK currently solves one frequency per run.
        std::vector<Eigen::VectorXd> zTabs(input.SSP.size());
        std::vector<Eigen::VectorXi> zFromSz(input.SSP.size());
        std::vector<Eigen::VectorXi> zFromRz(input.SSP.size());
        std::vector<int> nzTabs(input.SSP.size(), 0);

        int LRecordLength = std::max(2 * nfreq, 32);
        for (size_t iprof = 0; iprof < input.SSP.size(); ++iprof)
        {
            const auto &ssp = input.SSP.at(iprof);
            const int firstAc = ssp.FirstAcoustic;
            const int lastAc = ssp.LastAcoustic;
            if (firstAc < 0 || lastAc < firstAc || lastAc >= ssp.NMedia)
            {
                throw std::runtime_error("Cannot export MOD file: invalid acoustic media range.");
            }

            const auto &eigen = output.eigen[iprof];
            if (eigen.ModeZ.size() > 0)
            {
                zTabs[iprof] = eigen.ModeZ;
                nzTabs[iprof] = static_cast<int>(eigen.ModeZ.size());
                zFromSz[iprof].resize(0);
                zFromRz[iprof].resize(0);
            }
            else
            {
                MergeVectors(input.Pos.Sz, input.Pos.Rz, zTabs[iprof], nzTabs[iprof], zFromSz[iprof], zFromRz[iprof]);
            }
            const int nMedia = lastAc - firstAc + 1;
            LRecordLength = std::max(LRecordLength, 2 * nzTabs[iprof]);
            LRecordLength = std::max(LRecordLength, 3 * nMedia);
        }

        const int recordBytes = 4 * LRecordLength;
        std::ofstream MODFile(filename + ".mod", std::ios::binary | std::ios::trunc);
        if (!MODFile.is_open())
        {
            throw std::runtime_error("Failed to open MOD file for writing: " + filename + ".mod");
        }

        auto append_bytes = [](std::vector<char> &record, const void *data, size_t size)
        {
            const char *ptr = static_cast<const char *>(data);
            record.insert(record.end(), ptr, ptr + size);
        };
        auto append_int = [&](std::vector<char> &record, int value)
        {
            append_bytes(record, &value, sizeof(value));
        };
        auto append_float = [&](std::vector<char> &record, float value)
        {
            append_bytes(record, &value, sizeof(value));
        };
        auto append_double = [&](std::vector<char> &record, double value)
        {
            append_bytes(record, &value, sizeof(value));
        };
        auto append_fixed_string = [&](std::vector<char> &record, const std::string &value, size_t width)
        {
            std::string padded(width, ' ');
            std::memcpy(padded.data(), value.data(), std::min(width, value.size()));
            append_bytes(record, padded.data(), padded.size());
        };
        auto append_complex_float = [&](std::vector<char> &record, std::complex<double> value)
        {
            const float realPart = static_cast<float>(std::real(value));
            const float imagPart = static_cast<float>(std::imag(value));
            append_float(record, realPart);
            append_float(record, imagPart);
        };
        auto write_record = [&](int recordNumber, const std::vector<char> &payload)
        {
            if (payload.size() > static_cast<size_t>(recordBytes))
            {
                throw std::runtime_error("MOD record payload exceeds fixed record length.");
            }
            std::vector<char> record(recordBytes, 0);
            if (!payload.empty())
            {
                std::memcpy(record.data(), payload.data(), payload.size());
            }
            const std::streamoff offset = static_cast<std::streamoff>(recordNumber - 1) * recordBytes;
            MODFile.seekp(offset, std::ios::beg);
            MODFile.write(record.data(), static_cast<std::streamsize>(record.size()));
            if (!MODFile)
            {
                throw std::runtime_error("Failed while writing MOD file record.");
            }
        };
        auto bc_char = [](BC_Mode mode) -> char
        {
            switch (mode)
            {
            case BC_Mode::MODE_R_Rigid:
                return 'R';
            case BC_Mode::MODE_V_Vacuum:
                return 'V';
            case BC_Mode::MODE_F_File:
                return 'F';
            case BC_Mode::MODE_A_Half_space:
                return 'A';
            case BC_Mode::MODE_G_Grain:
                return 'G';
            case BC_Mode::MODE_P_Precomputed:
                return 'P';
            default:
                return 'V';
            }
        };
        auto material_name = [](Media_Mode mode) -> std::string
        {
            return mode == Media_Mode::MODE_E_Elastic ? "ELASTIC" : "ACOUSTIC";
        };
        auto hs_cp = [](const HSInfo &hs) -> std::complex<double>
        {
            return hs.cp != std::complex<double>(0.0, 0.0) ? hs.cp : std::complex<double>(hs.alphaR, hs.alphaI);
        };
        auto hs_cs = [](const HSInfo &hs) -> std::complex<double>
        {
            return hs.cs != std::complex<double>(0.0, 0.0) ? hs.cs : std::complex<double>(hs.betaR, hs.betaI);
        };
        auto medium_top_depth = [](const ssp::SSPStructure &ssp, int medium) -> float
        {
            if (ssp.offset.size() > medium && ssp.z.size() > ssp.offset(medium))
            {
                return static_cast<float>(ssp.z(ssp.offset(medium)));
            }
            return 0.0f;
        };
        auto medium_bottom_depth = [&](const ssp::SSPStructure &ssp, int medium) -> float
        {
            const int end = ssp.get_media_end(medium);
            if (ssp.z.size() > end)
            {
                return static_cast<float>(ssp.z(end));
            }
            return medium_top_depth(ssp, medium) + (ssp.depth.size() > medium ? static_cast<float>(ssp.depth(medium)) : 0.0f);
        };

        int iRecProfile = 1;
        for (size_t iprof = 0; iprof < input.SSP.size(); ++iprof)
        {
            const auto &ssp = input.SSP.at(iprof);
            const auto &trid = this->intm_TridMtx[iprof];
            const auto &eigen = output.eigen[iprof];
            const int firstAc = ssp.FirstAcoustic;
            const int lastAc = ssp.LastAcoustic;
            const int nMedia = lastAc - firstAc + 1;
            const int nzTab = nzTabs[iprof];
            const int m = std::max(0, eigen.M);

            std::vector<char> record;
            append_int(record, LRecordLength);
            append_fixed_string(record, input.Title, 80);
            append_int(record, nfreq);
            append_int(record, nMedia);
            append_int(record, nzTab);
            append_int(record, nzTab);
            write_record(iRecProfile, record);

            record.clear();
            for (int medium = firstAc; medium <= lastAc; ++medium)
            {
                const int nMesh = (trid.N.size() > medium && trid.N(medium) > 0) ? trid.N(medium) : ssp.NMesh(medium);
                append_int(record, nMesh);
                append_fixed_string(record, material_name(ssp.Material.at(medium)), 8);
            }
            write_record(iRecProfile + 1, record);

            record.clear();
            for (int medium = firstAc; medium <= lastAc; ++medium)
            {
                append_float(record, medium_top_depth(ssp, medium));
                const int loc = (trid.Loc.size() > medium) ? trid.Loc(medium) : ssp.offset(medium);
                const float rhoTop = (trid.rho.size() > loc && trid.rho(loc) != 0.0)
                                         ? static_cast<float>(trid.rho(loc))
                                         : static_cast<float>(ssp.rho(ssp.offset(medium)));
                append_float(record, rhoTop);
            }
            write_record(iRecProfile + 2, record);

            record.clear();
            append_double(record, input.freqinfo.freq);
            write_record(iRecProfile + 3, record);

            record.clear();
            for (int iz = 0; iz < nzTab; ++iz)
            {
                append_float(record, static_cast<float>(zTabs[iprof](iz)));
            }
            write_record(iRecProfile + 4, record);

            iRecProfile += 5;

            record.clear();
            append_int(record, m);
            write_record(iRecProfile, record);

            record.clear();
            append_fixed_string(record, std::string(1, bc_char(ssp.HSTop.BC)), 1);
            append_complex_float(record, hs_cp(ssp.HSTop));
            append_complex_float(record, hs_cs(ssp.HSTop));
            append_float(record, static_cast<float>(ssp.HSTop.rho));
            append_float(record, medium_top_depth(ssp, firstAc));
            append_fixed_string(record, std::string(1, bc_char(ssp.HSBot.BC)), 1);
            append_complex_float(record, hs_cp(ssp.HSBot));
            append_complex_float(record, hs_cs(ssp.HSBot));
            append_float(record, static_cast<float>(ssp.HSBot.rho));
            append_float(record, medium_bottom_depth(ssp, lastAc));
            write_record(iRecProfile + 1, record);

            for (int mode = 0; mode < m; ++mode)
            {
                record.clear();
                for (int iz = 0; iz < nzTab; ++iz)
                {
                    std::complex<double> phi = 0.0;
                    if (eigen.PhiMode.rows() > mode && eigen.PhiMode.cols() == nzTab)
                    {
                        phi = eigen.PhiMode(mode, iz);
                    }
                    else
                    {
                        int sourceIndex = -1;
                        for (int isz = 0; isz < zFromSz[iprof].size(); ++isz)
                        {
                            if (zFromSz[iprof](isz) == iz)
                            {
                                sourceIndex = isz;
                                break;
                            }
                        }
                        if (sourceIndex >= 0)
                        {
                            phi = eigen.PsiS(mode, sourceIndex);
                        }
                        else
                        {
                            for (int irz = 0; irz < zFromRz[iprof].size(); ++irz)
                            {
                                if (zFromRz[iprof](irz) == iz)
                                {
                                    phi = eigen.PsiR(mode, irz);
                                    break;
                                }
                            }
                        }
                    }
                    append_complex_float(record, phi);
                }
                write_record(iRecProfile + 2 + mode, record);
            }

            int iFirst = 0;
            const int modesPerRecord = std::max(1, LRecordLength / 2);
            const int kRecordCount = (m > 0) ? (1 + (2 * m - 1) / LRecordLength) : 0;
            for (int irec = 0; irec < kRecordCount; ++irec)
            {
                record.clear();
                const int iLast = std::min(m, iFirst + modesPerRecord);
                for (int mode = iFirst; mode < iLast; ++mode)
                {
                    append_complex_float(record, eigen.k(mode));
                }
                write_record(iRecProfile + m + 2 + irec, record);
                iFirst = iLast;
            }

            iRecProfile += 3 + m + ((m > 0) ? ((2 * m - 1) / LRecordLength) : 0);
        }

        MODFile.close();
    }

    void Interface::export_shd(std::string filename, int dataType)
    {
        // std::cout << "开始导出SHD文件: " << filename << ", 数据类型: " << dataType << std::endl;
        int LRecl;
        auto &output = this->getOutput(); // 获取输出
        auto &input = this->getParams();  // 获取输入
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
        std::vector<double> freqvec(Nfreq, input.freqinfo.freq);
        std::vector<float> theta(Ntheta, 0.0), Sx(NSx, 0.0), Sy(NSy, 0.0);
        std::vector<float> Sz(input.Pos.NSz);
        for (int i = 0; i < input.Pos.NSz; i++)
        {
            Sz[i] = input.Pos.Sz[i];
        }
        std::vector<float> Rz(input.Pos.NRz);
        for (int i = 0; i < input.Pos.NRz; i++)
        {
            Rz[i] = input.Pos.Rz[i];
        }
        std::vector<float> Rr(input.Pos.NRr);
        for (int i = 0; i < input.Pos.NRr; i++)
        {
            Rr[i] = input.Pos.Rr[i];
        }

        std::string plottype;
        switch (input.Pos.GridType)
        {
        case Grid_Mode::MODE_I_Irregular:
            plottype = "irregular ";
            break;
        case Grid_Mode::MODE_R_Rectangular:
            plottype = "rectilin  ";
            break;
        }

        LRecl = std::max(std::max(std::max(std::max(std::max(41, 2 * input.freqinfo.Nfreq), Ntheta), input.Pos.NSz), input.Pos.NRz), 2 * input.Pos.NRr);

        SHDFile.write(reinterpret_cast<char *>(&LRecl), sizeof(int));
        SHDFile.write(input.Title.c_str(), input.Title.size());
        SHDFile.seekp(1 * 4 * LRecl, std::ios::beg);
        SHDFile.write(plottype.c_str(), 80);
        SHDFile.seekp(2 * 4 * LRecl, std::ios::beg);
        SHDFile.write(reinterpret_cast<char *>(&Nfreq), sizeof(int));
        SHDFile.write(reinterpret_cast<char *>(&Ntheta), sizeof(int));
        SHDFile.write(reinterpret_cast<char *>(&NSx), sizeof(int));
        SHDFile.write(reinterpret_cast<char *>(&NSy), sizeof(int));
        SHDFile.write(reinterpret_cast<char *>(&input.Pos.NSz), sizeof(int));
        SHDFile.write(reinterpret_cast<char *>(&input.Pos.NRz), sizeof(int));
        SHDFile.write(reinterpret_cast<char *>(&input.Pos.NRr), sizeof(int));
        float freq = static_cast<float>(input.freqinfo.freq); // 转换double为float
        SHDFile.write(reinterpret_cast<char *>(&freq), sizeof(float));
        SHDFile.write(reinterpret_cast<char *>(&Atten), sizeof(float));
        SHDFile.seekp(3 * 4 * LRecl, std::ios::beg);
        SHDFile.write(reinterpret_cast<char *>(freqvec.data()), Nfreq * sizeof(double));
        SHDFile.seekp(4 * 4 * LRecl, std::ios::beg);
        SHDFile.write(reinterpret_cast<char *>(theta.data()), Ntheta * sizeof(float));
        SHDFile.seekp(5 * 4 * LRecl, std::ios::beg);
        SHDFile.write(reinterpret_cast<char *>(Sx.data()), 1 * sizeof(float));
        SHDFile.seekp(6 * 4 * LRecl, std::ios::beg);
        SHDFile.write(reinterpret_cast<char *>(Sy.data()), 1 * sizeof(float));
        SHDFile.seekp(7 * 4 * LRecl, std::ios::beg);
        SHDFile.write(reinterpret_cast<char *>(Sz.data()), input.Pos.NSz * sizeof(float));
        SHDFile.seekp(8 * 4 * LRecl, std::ios::beg);
        SHDFile.write((char *)Rz.data(), Rz.size() * sizeof(float));
        SHDFile.seekp(9 * 4 * LRecl, std::ios::beg);
        SHDFile.write(reinterpret_cast<char *>(Rr.data()), Rr.size() * sizeof(float));

        // 选择数据
        std::complex<float> *data = output.u_AllSources;
        switch (dataType)
        {
        case 1:
            data = output.u_AllSources;
            break;
        case 2:
            data = output.v_AllSources;
            break;
        case 3:
            data = output.h_AllSources;
            break;
        default:
            data = output.u_AllSources;
            break;
        }

        // std::cout << "开始写入数据，NSz=" << input.Pos->NSz << ", NRz_per_range=" << input.Pos->NRz_per_range << ", NRr=" << input.Pos->NRr << std::endl;
        // 遍历所有位置并写入声压值的实部和虚部
        for (int isz = 0; isz < input.Pos.NSz; ++isz)
        {
            for (int irz = 0; irz < input.Pos.NRz_per_range; ++irz)
            {
                int recnum = 10 + isz * input.Pos.NRz_per_range + irz;
                SHDFile.seekp(recnum * 4 * LRecl, std::ios::beg);
                for (int ir = 0; ir < input.Pos.NRr; ++ir)
                {
                    std::complex<float> &P = data[GetFieldAddr(isz, irz, ir, &input.Pos)];
                    SHDFile.write(reinterpret_cast<const char *>(&P), sizeof(P));
                }
            }
            // // 添加进度信息
            // if (isz % 10 == 0 || isz == input.Pos->NSz - 1) {
            //     std::cout << "进度: " << (isz + 1) << "/" << input.Pos->NSz << std::endl;
            // }
        }

        SHDFile.close();
        // std::cout << "SHD文件导出完成: " << filename << std::endl;
    }

    OOK_parameters &Interface::getParams() const // 获取参数的引用
    {
        if (!this->params)
        {
            throw std::runtime_error("OpenOceanKraken::Interface parameters have been freed.");
        }
        return *this->params;
    }
    const OOK_parameters &Interface::getParams_const() const // 获取参数的副本
    {
        if (!this->params)
        {
            throw std::runtime_error("OpenOceanKraken::Interface parameters have been freed.");
        }
        return *this->params;
    }
    OOK_output &Interface::getOutput() const // 获取输出的引用
    {
        if (!this->output)
        {
            throw std::runtime_error("OpenOceanKraken::Interface output has been freed.");
        }
        return *this->output;
    }
    OOK_output Interface::getOutput_Copy() const // 获取输出的副本
    {
        throw std::runtime_error("getOutput_Copy() is disabled because OOK_output owns raw buffers. Use getOutput_const() or typed getters instead.");
    }
    const OOK_output &Interface::getOutput_const() const // 获取输出的副本
    {
        if (!this->output)
        {
            throw std::runtime_error("OpenOceanKraken::Interface output has been freed.");
        }
        return *this->output;
    }

    void Interface::set_RProf(const Eigen::VectorXd &ranges)
    {
        ensureAlive();
        if (ranges.size() <= 0)
        {
            throw std::invalid_argument("Profile ranges must not be empty.");
        }
        for (Eigen::Index index = 0; index < ranges.size(); ++index)
        {
            if (!std::isfinite(ranges(index)))
            {
                throw std::invalid_argument("Profile ranges must be finite.");
            }
            if (index > 0 && ranges(index) < ranges(index - 1))
            {
                throw std::invalid_argument("Profile ranges must be non-decreasing.");
            }
        }
        auto &input = getParams();
        input.RProf = ranges;
        input.NProf = static_cast<int>(ranges.size());
        markDirty(DirtyKind::All);
    }

    void Interface::set_RProf(double start, double end, int count)
    {
        if (!std::isfinite(start) || !std::isfinite(end) || end < start)
        {
            throw std::invalid_argument("Profile range bounds must be finite and non-decreasing.");
        }
        if (count <= 0)
        {
            throw std::invalid_argument("Profile range count must be greater than zero.");
        }
        set_RProf(Eigen::VectorXd::LinSpaced(count, start, end));
    }

    void Interface::set_MLimit(int limit)
    {
        ensureAlive();
        if (limit <= 0)
        {
            throw std::invalid_argument("Mode limit must be greater than zero.");
        }
        getParams().MLimit = limit;
        markDirty(DirtyKind::All);
    }

    void Interface::set_CoherenceType(CoherenceType type)
    {
        ensureAlive();
        getParams().coherenceType = type;
        markDirty(DirtyKind::Field);
    }

    void Interface::set_ModeType(ModeType type)
    {
        ensureAlive();
        getParams().modeType = type;
        markDirty(DirtyKind::Field);
    }

    namespace
    {
        FieldSnapshot copy_field_snapshot(
            const OOK_parameters &input,
            const std::complex<float> *data,
            const char *label)
        {
            if (!data)
            {
                throw std::runtime_error(std::string(label) + " output is not available.");
            }
            FieldSnapshot snapshot;
            snapshot.title = input.Title;
            snapshot.frequency = input.freqinfo.freq;
            snapshot.grid_type = input.Pos.GridType;
            snapshot.source_count = static_cast<std::size_t>(input.Pos.NSz);
            snapshot.range_count = static_cast<std::size_t>(input.Pos.NRr);
            snapshot.depth_count = static_cast<std::size_t>(input.Pos.NRz_per_range);
            snapshot.source_depths.assign(input.Pos.Sz.data(), input.Pos.Sz.data() + input.Pos.Sz.size());
            snapshot.receiver_ranges.assign(input.Pos.Rr.data(), input.Pos.Rr.data() + input.Pos.Rr.size());
            snapshot.receiver_depths.assign(input.Pos.Rz.data(), input.Pos.Rz.data() + input.Pos.Rz.size());
            const std::size_t value_count = snapshot.source_count * snapshot.range_count * snapshot.depth_count;
            snapshot.values.assign(data, data + value_count);
            return snapshot;
        }
    }

    FieldSnapshot Interface::getPressureCopy() const
    {
        ensureSetup();
        return copy_field_snapshot(getParams_const(), getOutput_const().u_AllSources, "Pressure");
    }

    FieldSnapshot Interface::getVerticalVelocityCopy() const
    {
        ensureSetup();
        return copy_field_snapshot(getParams_const(), getOutput_const().v_AllSources, "Vertical velocity");
    }

    FieldSnapshot Interface::getHorizontalVelocityCopy() const
    {
        ensureSetup();
        return copy_field_snapshot(getParams_const(), getOutput_const().h_AllSources, "Horizontal velocity");
    }

    std::vector<ModeProfileSnapshot> Interface::getModesCopy() const
    {
        ensureAlive();
        const auto &input = getParams_const();
        const auto &result = getOutput_const();
        if (!result.eigen)
        {
            throw std::runtime_error("Eigenmode output is not available. Call runEigen() or run() first.");
        }

        std::vector<ModeProfileSnapshot> snapshots;
        snapshots.reserve(static_cast<std::size_t>(input.NProf));
        for (int profile = 0; profile < input.NProf; ++profile)
        {
            const auto &eigen = result.eigen[profile];
            if (eigen.M < 0 || eigen.k.size() < eigen.M)
            {
                throw std::runtime_error("Eigenmode output is incomplete for profile " + std::to_string(profile) + ".");
            }

            ModeProfileSnapshot snapshot;
            if (profile < input.RProf.size())
            {
                snapshot.profile_range = input.RProf(profile);
            }
            snapshot.wavenumbers.assign(eigen.k.data(), eigen.k.data() + eigen.M);
            if (eigen.VG.size() >= eigen.M)
            {
                snapshot.group_velocity.assign(eigen.VG.data(), eigen.VG.data() + eigen.M);
            }

            if (eigen.ModeZ.size() > 0 &&
                eigen.PhiMode.rows() >= eigen.M &&
                eigen.PhiMode.cols() == eigen.ModeZ.size())
            {
                snapshot.depth.assign(eigen.ModeZ.data(), eigen.ModeZ.data() + eigen.ModeZ.size());
                snapshot.mode_shapes = eigen.PhiMode.topRows(eigen.M);
            }
            else if (eigen.PsiR.rows() >= eigen.M && eigen.PsiR.cols() > 0)
            {
                const Eigen::Index depth_count = std::min(eigen.PsiR.cols(), input.Pos.Rz.size());
                snapshot.depth.assign(input.Pos.Rz.data(), input.Pos.Rz.data() + depth_count);
                snapshot.mode_shapes = eigen.PsiR.topLeftCorner(eigen.M, depth_count);
            }
            else
            {
                snapshot.mode_shapes.resize(eigen.M, 0);
            }
            snapshots.push_back(std::move(snapshot));
        }
        return snapshots;
    }

}
