
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
          impl(std::make_unique<OpenOceanKraken_PIMPL>())
    {
        this->init(); // 初始化
    }
    Interface::Interface(ThreadPool &pool)
        : params(new OOK_parameters()),
          output(new OOK_output()),
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
        std::cout << "OpenOcean-Kraken: Initialized" << std::endl;
    }

    void Interface::input_setup()
    {
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
        auto &params = this->getParams(); // 获取参数
        size_t num_threads = this->NumThreads;
        this->intm_TridMtx = new TridMtx[num_threads];
        for (size_t i = 0; i < num_threads; i++)
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
        this->NumThreads = num_threads;
    }
    void Interface::setThreadPool(ThreadPool &pool)
    {
        this->threadPool = &pool;
    }

    int Interface::getNumThreads() const
    {
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
        auto &params = this->getParams(); // 获取参数
        auto &output = this->getOutput(); // 获取输出
        for (size_t iprof = 0; iprof < params.SSP.size(); iprof++)
        {
            FieldWorker(iprof, params, output);
        }
    }

    void Interface::runEigen() // 运行特征值
    {
        auto &params = this->getParams(); // 获取参数
        auto &output = this->getOutput(); // 获取输出
        for (size_t iprof = 0; iprof < params.SSP.size(); iprof++)
        {
            EigenVWorker(*this->threadPool, this->NumThreads, iprof, params, *this->intm_TridMtx, output);
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
        this->impl->OUTPUT_FIELD.Finalize(*this->output);
        this->impl->OUTPUT_EIGEN.Finalize(*this->output);
        delete params;
        delete this->output;
        if (this->intm_TridMtx != nullptr)
        {
            delete[] this->intm_TridMtx;
        }

        std::cout << "free memory successfully!" << std::endl;
    }

    // 参数设置
    void Interface::set_Title(std::string &title)
    {
        auto &params = this->getParams(); // 获取参数
        this->impl->INPUT_FREQ.set_title(params, title);
    }
    void Interface::set_Freq(double freq)
    {
        auto &params = this->getParams(); // 获取参数
        this->impl->INPUT_FREQ.set_freq(params, freq);
    }
    void Interface::set_freqvec(Eigen::VectorXd freqvec) // 设置频率向量
    {
        auto &params = this->getParams(); // 获取参数
        this->impl->INPUT_FREQ.set_freqvec(params, freqvec);
    }
    void Interface::set_SSP(const std::vector<ssp::Range_Independent_Area> &sspInput) // 设置SSP
    {
        auto &params = this->getParams(); // 获取参数
        this->impl->INPUT_SSP.set_SSP(params, sspInput);
    }

    void Interface::set_AttenUnit(Atten_Mode mode) // 设置衰减单位
    {
        auto &params = this->getParams(); // 获取参数
        this->impl->INPUT_SSP.set_AttenUnit(params, mode);
    }

    void Interface::set_Sz(const Eigen::VectorXd &Sz)
    {
        auto &params = this->getParams(); // 获取参数
        this->impl->INPUT_SZ_RZ_RR.set_Sz(params, Sz);
    }
    void Interface::set_Sz(const double &start, const double &end, const int &NSz) // 设置声源深度（插值）
    {
        auto &params = this->getParams(); // 获取参数
        this->impl->INPUT_SZ_RZ_RR.set_Sz(params, start, end, NSz);
    }
    void Interface::set_Rr(const Eigen::VectorXd &Rr) // 设置水平接收
    {
        auto &params = this->getParams(); // 获取参数
        this->impl->INPUT_SZ_RZ_RR.set_Rr(params, Rr);
    }
    void Interface::set_Rr(const double &start, const double &end, const int &NRr) // 设置水平接收（插值）
    {
        auto &params = this->getParams(); // 获取参数
        this->impl->INPUT_SZ_RZ_RR.set_Rr(params, start, end, NRr);
    }
    void Interface::set_Rz(const Eigen::VectorXd &Rz)
    {
        auto &params = this->getParams(); // 获取参数
        this->impl->INPUT_SZ_RZ_RR.set_Rz(params, Rz);
    }
    void Interface::set_Rz(const double &start, const double &end, const int &NRz) // 设置接收深度（插值）
    {
        auto &params = this->getParams(); // 获取参数
        this->impl->INPUT_SZ_RZ_RR.set_Rz(params, start, end, NRz);
    }
    void Interface::set_Ro(const Eigen::VectorXd &Ro) // 设置阵列倾斜
    {
        auto &params = this->getParams(); // 获取参数
        this->impl->INPUT_SZ_RZ_RR.set_Ro(params, Ro);
    }
    void Interface::set_Ro(const double &start, const double &end, const int &NRo) // 设置阵列倾斜
    {
        auto &params = this->getParams(); // 获取参数
        this->impl->INPUT_SZ_RZ_RR.set_Ro(params, start, end, NRo);
    }

    void Interface::set_cPhase(double cLow, double cHigh) // 设置最低频率
    {
        auto &params = this->getParams(); // 获取参数
        this->impl->INPUT_FREQ.set_cPhase(params, cLow, cHigh);
    }
    void Interface::set_GridType(Grid_Mode type) // 设置网格类型
    {
        auto &params = this->getParams(); // 获取参数
        this->impl->INPUT_SZ_RZ_RR.set_GridType(params, type);
    }

    void Interface::set_Rmax(double Rmax) // 设置最大计算距离，用于缩放error
    {
        auto &params = this->getParams(); // 获取参数
        this->impl->INPUT_SZ_RZ_RR.set_RMax(params, Rmax);
    }
    void Interface::set_SourceType(Source_Mode type) // 设置源类型
    {
        auto &params = this->getParams(); // 获取参数
        this->impl->INPUT_FREQ.set_SourceType(params, type);
    }
    void Interface::set_RunMode(Run_Mode mode)
    {
        auto &params = this->getParams(); // 获取参数
        this->impl->INPUT_FREQ.set_RunMode(params, mode);
    } // 运行模式
    void Interface::set_Velocity_enable(bool is_Velocity) // 设置是否计算振速
    {
        auto &params = this->getParams(); // 获取参数
        this->impl->INPUT_FREQ.set_Velocity_enable(params, is_Velocity);
    }
    void Interface::set_ReflCoef_Top(std::vector<ReflectionCoef> ReflCoef)
    {
        auto &params = this->getParams(); // 获取参数
        this->impl->INPUT_REFLCOEF.set_ReflCoef_Top(params, ReflCoef);
    } // 设置顶部反射系数
    void Interface::set_ReflCoef_Bottom(std::vector<ReflectionCoef> ReflCoef)
    {
        auto &params = this->getParams(); // 获取参数
        this->impl->INPUT_REFLCOEF.set_ReflCoef_Bot(params, ReflCoef);
    } // 设置底部反射系数
    void Interface::set_SBP(const Eigen::VectorXd &pat, const Eigen::VectorXd &theta) // 设置指向性
    {
        auto &params = this->getParams(); // 获取参数
        this->impl->INPUT_SBP.set_Pat(params, pat, theta);
    }

    bool Interface::to_json(const std::string &jsonPath) const // 将参数写入json
    {


    }
    bool Interface::from_json(const std::string &jsonPath) // 从json读取参数
    {


    }


    std::string Interface::to_json_string() const // 将参数写入json字符串
    {


    }

    std::complex<float> *Interface::get_u(int srcIndex) // 获取某个声源复声压指针
    {
        auto &output = this->getOutput_const(); // 获取输出
        auto &input = this->getParams_const();  // 获取输入
        return output.u_AllSources + srcIndex * input.Pos.NRr * input.Pos.NRz_per_range;
    }

    std::complex<float> *Interface::get_v(int srcIndex) // 获取某个声源复垂直振速
    {
        auto &output = this->getOutput_const(); // 获取输出
        auto &input = this->getParams_const();  // 获取输入
        return output.v_AllSources + srcIndex * input.Pos.NRr * input.Pos.NRz_per_range;
    }

    std::complex<float> *Interface::get_h(int srcIndex) // 获取某个声源水平振速
    {
        auto &output = this->getOutput_const(); // 获取输出
        auto &input = this->getParams_const();  // 获取输入
        return output.h_AllSources + srcIndex * input.Pos.NRr * input.Pos.NRz_per_range;
    }

    std::complex<float> *Interface::get_u_AllSources() // 获取全部声源的复声压
    {
        auto &output = this->getOutput_const(); // 获取输出
        return output.u_AllSources;
    }

    std::complex<float> *Interface::get_v_AllSources() // 获取全部声源的垂直振速
    {
        auto &output = this->getOutput_const(); // 获取输出
        return output.v_AllSources;
    }

    std::complex<float> *Interface::get_h_AllSources() // 获取全部声源的水平振速
    {
        auto &output = this->getOutput_const(); // 获取输出
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

    OOK_parameters &Interface::getParams() // 获取参数的引用
    {
        return *this->params;
    }
    const OOK_parameters &Interface::getParams_const() const // 获取参数的副本
    {
        return *this->params;
    }
    OOK_output &Interface::getOutput() const // 获取输出的引用
    {
        return *this->output;
    }
    OOK_output Interface::getOutput_Copy() const // 获取输出的副本
    {
        return *this->output;
    }
    const OOK_output &Interface::getOutput_const() const // 获取输出的副本
    {
        return *this->output;
    }

}
