#include "kkc_log.h"
#include "kkc_interface.h"
#include "util.h"
#include "input_boundary.hpp"
#include "input_SSP.hpp"
#include "input_Sz_Rz_RR.hpp"
#include "input_reflcoef.hpp"
#include "input_sbp.hpp"
#include "input_Freq.hpp"
#include "output_eigen.hpp"
#include "output_field.hpp"
#include "run.h"

void print_header()
{
    std::cout << " - OpenOcean-Kraken: A parallel underwater acoustic nomal modes simulator with pybind11 Python and C++ interfaces\n"
                 "\n"
                 "Copyright (C) 2024-2025 OpenOceanAcoustic\n"
                 "Authors: Qian Peng, Dai Nuoge, Liang Yi, Luo Shixiong, Liu Ruihang, Liu Jiatong, Wang Zhengwei\n"
                 "Mail: yingxinliang1@gmail.com\n"
                 "Based on Kraken, which is Copyright (C) 1983-2024 Michael B. Porter\n"
                 "GPL3 licensed, no warranty, see LICENSE or https://www.gnu.org/licenses/\n"
                 "\n";
}

// 定义实现类
class kkc_interface_PIMPL
{
public:
    input_Freq INPUT_FREQ;
    input_SSP INPUT_SSP;
    input_Boundary INPUT_BOUNDARY;
    input_Sz_Rz_RR INPUT_SZ_RZ_RR;
    input_reflcoef INPUT_REFLCOEF;
    input_sbp INPUT_SBP;

    output_Eigen OUTPUT_EIGEN;
    output_Field OUTPUT_FIELD;
    kkc_interface_PIMPL()
    {
    }
    ~kkc_interface_PIMPL()
    {
    }
};

// 默认构造函数定义，不接收线程池参数
kkc_interface::kkc_interface()
    : params(new parameters()),
      output(new kkc_output()),
      impl(std::make_unique<kkc_interface_PIMPL>())
{
    this->init(); // 初始化
}

kkc_interface::~kkc_interface()
{
    // 释放params的各个成员变量的内存
    this->free(); // 释放内存
}

void kkc_interface::init() // 初始化
{
    if (!params)
    {
        std::cerr << "Error: Failed to allocate memory for parameters" << std::endl;
        exit(1);
    }

    // 初始化params的各个成员变量
    params->NProf = 1;
    params->freqinfo = std::make_unique<FreqInfo>();
    params->SSP = std::make_unique<SSPStructure[]>(params->NProf);
    params->HSTop = std::make_unique<HSInfo[]>(params->NProf);
    params->HSBot = std::make_unique<HSInfo[]>(params->NProf);
    params->Pos = std::make_unique<Position>();
    params->SBP = std::make_unique<SrcBmPat>();
    params->log = &kkc_Log::get_instance();
    this->field_size = 0; // 场大小初始化为0

    this->is_setup = false;           // 初始化是否设置参数标志位
    auto &params = this->getParams(); // 获取参数
    auto &output = this->getOutput(); // 获取输出

    // 参数初始化
    impl->INPUT_FREQ.Init(params);
    impl->INPUT_SSP.Init(params);
    impl->INPUT_BOUNDARY.Init(params);
    impl->INPUT_SZ_RZ_RR.Init(params);
    impl->INPUT_REFLCOEF.Init(params);
    impl->INPUT_SBP.Init(params);

    // default params
    impl->INPUT_FREQ.Default(params);
    impl->INPUT_SSP.Default(params);
    impl->INPUT_BOUNDARY.Default(params);
    impl->INPUT_SZ_RZ_RR.Default(params);
    impl->INPUT_REFLCOEF.Default(params);
    impl->INPUT_SBP.Default(params);

    // output初始化
    impl->OUTPUT_FIELD.Init(output);
}

// 配置
void kkc_interface::setup() // 设置参数
{
    input_setup(); // 设置输入参数(外界输入参数，因此每次计算时候都重新计算一次环境)
    intm_setup();  // 设置中间矩阵参数

    if (!this->is_setup)
    {
        this->output_setup();  // 若未配置，则配置输出的内存
        this->is_setup = true; // 设置配置标志位
    }
}

// 输入配置
void kkc_interface::input_setup() // 设置输入参数
{
    auto &params = this->getParams(); // 获取参数
    auto &output = this->getOutput(); // 获取输出
    // 预处理
    impl->INPUT_SSP.Preprocess(params);
    impl->INPUT_BOUNDARY.Preprocess(params);
    impl->INPUT_SZ_RZ_RR.Preprocess(params);
    impl->INPUT_REFLCOEF.Preprocess(params);
    impl->INPUT_SBP.Preprocess(params);
    impl->INPUT_FREQ.Preprocess(params);

    // 初始化声压共享内存区域
    this->field_size = (size_t)params.Pos->NSz * (size_t)params.Pos->NRz_per_range * (size_t)params.Pos->NRr;
}

// 中间矩阵配置
void kkc_interface::intm_setup() // 设置中间矩阵参数
{
    auto &params = this->getParams(); // 获取参数
    size_t num_threads = 1;
    this->intm_TridMtx = new TridMtx[num_threads];
    for (size_t i = 0; i < num_threads; i++)
    {
        this->intm_TridMtx[i].resize(params.NMeshMax, params.NMediaMax, params.mesh.NSets);
    }
}

// 输出配置
void kkc_interface::output_setup() // 设置输出参数
{
    auto &params = this->getParams(); // 获取参数
    auto &output = this->getOutput(); // 获取输出
    impl->OUTPUT_EIGEN.Preprocess(params, output);
    impl->OUTPUT_FIELD.Preprocess(params, output);
}

void kkc_interface::run() // 运行
{
    print_header();    // 打印头信息
    this->setup();     // 配置
    this->runSolveV(); // 求解本征值和本征函数
    this->runField();  // 求解声场
}

void kkc_interface::clearResults() // 清除结果
{
    auto &output = this->getOutput(); // 获取输出
    auto &params = this->getParams(); // 获取参数

    this->impl->OUTPUT_FIELD.ClearResults(params, output);
}

void kkc_interface::free() // 释放内存
{
    // 释放params的各个成员变量的内存
    delete params;

    // 释放output的各个成员变量的内存
    this->impl->OUTPUT_FIELD.Finalize(*this->output);
    delete this->output;
    std::cout << "free memory successfully!" << std::endl;
}

void kkc_interface::runSolveV()
{
    auto &paramsRef = this->getParams_const(); // 获取参数
    auto &output = this->getOutput();          // 获取输出

    for (int iprof = 0; iprof < paramsRef.NProf; iprof++)
    {
        EigenVWorker(iprof, paramsRef, this->intm_TridMtx[0], output);
    }
}

void kkc_interface::runField()
{
    // FieldSolveWorker();
    auto &params = this->getParams(); // 获取参数
    auto &output = this->getOutput(); // 获取输出
    for (int iprof = 0; iprof < params.NProf; iprof++)
    {
        FieldWorker(iprof, params, output);
    }
}

void kkc_interface::set_Title(std::string &title) // 设置标题
{
    auto &params = this->getParams();          // 获取参数
    impl->INPUT_FREQ.set_title(params, title); // 设置标题
}

void kkc_interface::set_Freq(double freq) // 设置频率
{
    auto &params = this->getParams();        // 获取参数
    impl->INPUT_FREQ.set_freq(params, freq); // 设置频率
}

void kkc_interface::set_freqvec(VectorXd freqvec) // 设置频率向量
{
    auto &params = this->getParams();              // 获取参数
    impl->INPUT_FREQ.set_freqvec(params, freqvec); // 设置频率向量
}

void kkc_interface::set_RProf(const VectorXd &RProf) // 设置距离剖面
{
    auto &params = this->getParams();          // 获取参数
    impl->INPUT_FREQ.set_RProf(params, RProf); // 设置 RProf
}

void kkc_interface::set_RProf(const double &start, const double &end, const int &NProf) // 设置距离剖面（插值）
{
    auto &params = this->getParams();                      // 获取参数
    impl->INPUT_FREQ.set_RProf(params, start, end, NProf); // 设置 RProf
}


//传vector
void kkc_interface::set_SSP(const std::vector<SSP_1D> &sspInput) // 设置SSP
{
    auto &params = this->getParams();                 // 获取参数
    impl->INPUT_SSP.set_SSP(params, sspInput); // 设置1D SSP
}


void kkc_interface::set_AttenUnit(Atten_Mode mode) // 设置衰减单位
{
    auto &params = this->getParams();            // 获取参数
    impl->INPUT_SSP.set_AttenUnit(params, mode); // 设置衰减单位
}

void kkc_interface::set_Sz(const VectorXd &Sz)
{
    auto &params = this->getParams(); // 获取参数
    auto &output = this->getOutput(); // 获取输出
    this->impl->OUTPUT_FIELD.Finalize(output);
    this->is_setup = false;
    impl->INPUT_SZ_RZ_RR.set_Sz(params, Sz); // 设置 Sz
}

void kkc_interface::set_Sz(const double &start, const double &end, const int &NSz) // 设置声源深度（插值）
{
    auto &params = this->getParams(); // 获取参数
    auto &output = this->getOutput(); // 获取输出
    this->impl->OUTPUT_FIELD.Finalize(output);
    this->is_setup = false;
    impl->INPUT_SZ_RZ_RR.set_Sz(params, start, end, NSz); // 设置 Sz
}

void kkc_interface::set_Rr(const VectorXd &Rr) // 设置水平接收
{
    auto &params = this->getParams(); // 获取参数
    auto &output = this->getOutput(); // 获取输出
    this->impl->OUTPUT_FIELD.Finalize(output);
    this->is_setup = false;
    impl->INPUT_SZ_RZ_RR.set_Rr(params, Rr); // 设置 Rr
}

void kkc_interface::set_Rr(const double &start, const double &end, const int &NRr) // 设置水平接收（插值）
{
    auto &params = this->getParams(); // 获取参数
    auto &output = this->getOutput(); // 获取输出
    this->impl->OUTPUT_FIELD.Finalize(output);
    this->is_setup = false;
    impl->INPUT_SZ_RZ_RR.set_Rr(params, start, end, NRr); // 设置 Rr
}

void kkc_interface::set_Rz(const VectorXd &Rz) // 设置垂直接收
{
    auto &params = this->getParams(); // 获取参数
    auto &output = this->getOutput(); // 获取输出
    this->impl->OUTPUT_FIELD.Finalize(output);
    this->is_setup = false;
    impl->INPUT_SZ_RZ_RR.set_Rz(params, Rz); // 设置 Rz
}

void kkc_interface::set_Rz(const double &start, const double &end, const int &NRz) // 设置垂直接收（插值）
{
    auto &params = this->getParams(); // 获取参数
    auto &output = this->getOutput(); // 获取输出
    this->impl->OUTPUT_FIELD.Finalize(output);
    this->is_setup = false;
    impl->INPUT_SZ_RZ_RR.set_Rz(params, start, end, NRz); // 设置 Rz
}

// 边界条件设置类
void kkc_interface::set_surface_Type(BC_Mode bc, size_t iprof) // 设置边界条件类型
{
    auto &params = this->getParams();                         // 获取参数
    impl->INPUT_BOUNDARY.set_surface_Type(params, bc, iprof); // 设置边界条件类型
}

void kkc_interface::set_bottom_Type(BC_Mode bc, size_t iprof) // 设置底部边界条件类型
{
    auto &params = this->getParams();                        // 获取参数
    impl->INPUT_BOUNDARY.set_bottom_Type(params, bc, iprof); // 设置底部边界条件类型
}

void kkc_interface::set_BottomLine(double zTemp, double alphaR, double alphaI, double betaR, double betaI, double rho, size_t iprof) // 设置底部半空间
{
    auto &params = this->getParams();                                                            // 获取参数
    impl->INPUT_BOUNDARY.setBottomLine(params, zTemp, alphaR, alphaI, betaR, betaI, rho, iprof); // 设置底部半空间
}

void kkc_interface::set_SurfaceLine(double zTemp, double alphaR, double alphaI, double betaR, double betaI, double rho, size_t iprof) // 设置表面半空间
{
    auto &params = this->getParams();                                                             // 获取参数
    impl->INPUT_BOUNDARY.setSurfaceLine(params, zTemp, alphaR, alphaI, betaR, betaI, rho, iprof); // 设置表面半空间
}

void kkc_interface::set_cPhase(double cLow, double cHigh) // 设置最低频率
{
    auto &params = this->getParams();                 // 获取参数
    impl->INPUT_FREQ.set_cPhase(params, cLow, cHigh); // 设置最低频率
}

void kkc_interface::set_GridType(Grid_Mode type) // 设置网格类型
{
    auto &params = this->getParams();                // 获取参数
    impl->INPUT_SZ_RZ_RR.set_GridType(params, type); // 设置网格类型
}

void kkc_interface::set_SourceType(Source_Mode type) // 设置源类型
{
    auto &params = this->getParams();              // 获取参数
    impl->INPUT_FREQ.set_SourceType(params, type); // 设置源类型
}

void kkc_interface::set_RunMode(Run_Mode mode) // 运行模式
{
    auto &params = this->getParams();           // 获取参数
    impl->INPUT_FREQ.set_RunMode(params, mode); // 设置运行模式
}

void kkc_interface::set_Velocity_enable(bool is_Velocity) // 设置是否计算振速
{
    auto &params = this->getParams(); // 获取参数
    params.is_Velocity = is_Velocity; // 设置是否计算振速
}

// 反射系数
void kkc_interface::set_ReflCoef_Top(std::vector<ReflectionCoef> ReflCoef) // 设置顶部反射系数
{
    auto &params = this->getParams();                        // 获取参数
    impl->INPUT_REFLCOEF.set_ReflCoef_Top(params, ReflCoef); // 设置顶部反射系数
}

void kkc_interface::set_ReflCoef_Bottom(std::vector<ReflectionCoef> ReflCoef) // 设置底部反射系数
{
    auto &params = this->getParams();                        // 获取参数
    impl->INPUT_REFLCOEF.set_ReflCoef_Bot(params, ReflCoef); // 设置底部反射系数
}

// 指向性
void kkc_interface::set_SBP(const VectorXd &pat, const VectorXd &theta) // 设置指向性
{
    auto &params = this->getParams();            // 获取参数
    impl->INPUT_SBP.set_Pat(params, pat, theta); // 设置指向性
}

VectorXd kkc_interface::get_BottomLine(size_t iprof)
{
    return impl->INPUT_BOUNDARY.get_BottomLine(this->getParams_const(), iprof);
}
VectorXd kkc_interface::get_SurfaceLine(size_t iprof)
{
    return impl->INPUT_BOUNDARY.get_SurfaceLine(this->getParams_const(), iprof);
}
double kkc_interface::get_freq()
{
    return impl->INPUT_FREQ.get_freq(this->getParams_const());
}
std::vector<ReflectionCoef> kkc_interface::get_ReflCoef_Top()
{
    return impl->INPUT_REFLCOEF.get_ReflCoef_Top(this->getParams_const());
}
std::vector<ReflectionCoef> kkc_interface::get_ReflCoef_Bottom()
{
    return impl->INPUT_REFLCOEF.get_ReflCoef_Bot(this->getParams_const());
}
std::pair<VectorXd, VectorXd> kkc_interface::get_SBP()
{
    return impl->INPUT_SBP.get_SBP(this->getParams_const());
}


VectorXd kkc_interface::get_Sz()
{
    return impl->INPUT_SZ_RZ_RR.get_Sz(this->getParams_const());
}
VectorXd kkc_interface::get_Rr()
{
    return impl->INPUT_SZ_RZ_RR.get_Rr(this->getParams_const());
}
VectorXd kkc_interface::get_Rz()
{
    return impl->INPUT_SZ_RZ_RR.get_Rz(this->getParams_const());
}

std::complex<float> *kkc_interface::get_u(int srcIndex) // 获取某个声源复声压指针
{
    auto &output = this->getOutput_const(); // 获取输出
    auto &input = this->getParams_const();  // 获取输入
    return output.u_AllSources + srcIndex * input.Pos->NRr * input.Pos->NRz_per_range;
}

std::complex<float> *kkc_interface::get_v(int srcIndex) // 获取某个声源复垂直振速
{
    auto &output = this->getOutput_const(); // 获取输出
    auto &input = this->getParams_const();  // 获取输入
    return output.v_AllSources + srcIndex * input.Pos->NRr * input.Pos->NRz_per_range;
}

std::complex<float> *kkc_interface::get_h(int srcIndex) // 获取某个声源水平振速
{
    auto &output = this->getOutput_const(); // 获取输出
    auto &input = this->getParams_const();  // 获取输入
    return output.h_AllSources + srcIndex * input.Pos->NRr * input.Pos->NRz_per_range;
}

std::complex<float> *kkc_interface::get_u_AllSources() // 获取全部声源的复声压
{
    auto &output = this->getOutput_const(); // 获取输出
    return output.u_AllSources;
}

std::complex<float> *kkc_interface::get_v_AllSources() // 获取全部声源的垂直振速
{
    auto &output = this->getOutput_const(); // 获取输出
    return output.v_AllSources;
}

std::complex<float> *kkc_interface::get_h_AllSources() // 获取全部声源的水平振速
{
    auto &output = this->getOutput_const(); // 获取输出
    return output.h_AllSources;
}

// 输出

void kkc_interface::export_result(std::string filename) // 导出结果到文件
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

void kkc_interface::export_mod(std::string filename) // 导出模型到文件
{
    auto &output = this->getOutput(); // 获取输出
    auto &params = this->getParams(); // 获取输入
    auto &eigen = output.eigen[0];
    int LRecl, iRecProf = 0;
    int NzTab = params.Pos->NRz;
    std::ofstream MODFile = std::ofstream(filename + ".mod", std::ios::binary);
    if (!MODFile.is_open())
    {
        std::cerr << "无法打开MOD文件: " << filename << ".mod" << std::endl;
        return;
    }

    int ifreq = 0, iprof = 0;
    // if (ifreq == 0 && iprof == 0)
    // {
    LRecl = std::max(2 * params.freqinfo->Nfreq, std::max(2 * NzTab, std::max(32, 3 * (params.SSP[0].LastAcoustic - params.SSP[0].FirstAcoustic + 1))));
    // }
    if (ifreq == 0)
    {
        iRecProf = 0;
        MODFile.seekp(iRecProf * 4 * LRecl, std::ios::beg);
        MODFile.write(reinterpret_cast<char *>(&LRecl), sizeof(int));
        // 固定title为80个char
        char title[80] = {0};
        strncpy(title, params.Title.c_str(), params.Title.size());
        MODFile.write(title, sizeof(title));
        MODFile.write(reinterpret_cast<char *>(&params.freqinfo->Nfreq), sizeof(int));
        int numAcoustic = params.SSP[0].LastAcoustic - params.SSP[0].FirstAcoustic + 1;
        MODFile.write(reinterpret_cast<char *>(&numAcoustic), sizeof(int));
        MODFile.write(reinterpret_cast<char *>(&NzTab), sizeof(int));
        MODFile.write(reinterpret_cast<char *>(&NzTab), sizeof(int));

        MODFile.seekp((iRecProf + 1) * 4 * LRecl, std::ios::beg);
        for (int im = params.SSP[iprof].FirstAcoustic; im <= params.SSP[iprof].LastAcoustic; im++)
        {
            MODFile.write(reinterpret_cast<char *>(&params.SSP[iprof].NMesh(im)), sizeof(int));
            MODFile.write(reinterpret_cast<char *>(&params.SSP[iprof].Material[im]), sizeof(params.SSP[iprof].Material[im])); // fortran是字符串，这里是整数。而且只有一层，应该是params.SSP->Material[im]?
        }

        MODFile.seekp((iRecProf + 2) * 4 * LRecl, std::ios::beg);
        for (int im = params.SSP[iprof].FirstAcoustic; im <= params.SSP[iprof].LastAcoustic; im++)
        {
            float depth = static_cast<float>(params.SSP[iprof].depth[im]); // 转换double为float
            MODFile.write(reinterpret_cast<char *>(&depth), sizeof(float));
            float rho = static_cast<float>(params.SSP[iprof].rho[params.SSP[iprof].get_media_start(im)]); // 转换double为float
            MODFile.write(reinterpret_cast<char *>(&rho), sizeof(float));
        }
        MODFile.seekp((iRecProf + 3) * 4 * LRecl, std::ios::beg);
        MODFile.write(reinterpret_cast<char *>(params.freqinfo->freqvec.data()), params.freqinfo->Nfreq * sizeof(double));
        MODFile.seekp((iRecProf + 4) * 4 * LRecl, std::ios::beg);
        // zTab转为float输出
        VectorXf zTabf = params.Pos->Rz.cast<float>();
        MODFile.write(reinterpret_cast<char *>(zTabf.data()), NzTab * sizeof(float));
        // for (int isz = 0; isz < params.Pos->NSz; isz++)
        // {
        //     params.MODFile.write(reinterpret_cast<char*>(&params.Pos->Sz(isz)), sizeof(double));
        // }
        // for (int irz = 0; irz < params.Pos->NRz; irz++)
        // {
        //     params.MODFile.write(reinterpret_cast<char*>(&params.Pos->Rz(irz)), sizeof(double));
        // }
        iRecProf += 5;
    }
    MODFile.seekp((iRecProf + 1) * 4 * LRecl, std::ios::beg);
    MODFile.write(reinterpret_cast<char *>(&params.HSTop[iprof].BC), sizeof(params.HSTop[iprof].BC)); // 这也是整数而不是字符串
    MODFile.write(reinterpret_cast<char *>(&params.HSTop[iprof].cp), sizeof(std::complex<double>));
    MODFile.write(reinterpret_cast<char *>(&params.HSTop[iprof].cs), sizeof(std::complex<double>));
    float HSToprho = static_cast<float>(params.HSTop[iprof].rho); // 转换double为float
    MODFile.write(reinterpret_cast<char *>(&HSToprho), sizeof(float));
    float SSPz0 = static_cast<float>(params.SSP[iprof].z(0)); // 转换double为float
    MODFile.write(reinterpret_cast<char *>(&SSPz0), sizeof(float));
    MODFile.write(reinterpret_cast<char *>(&params.HSBot[iprof].BC), sizeof(params.HSBot[iprof].BC)); // 这也是整数而不是字符串
    MODFile.write(reinterpret_cast<char *>(&params.HSBot[iprof].cp), sizeof(std::complex<double>));
    MODFile.write(reinterpret_cast<char *>(&params.HSBot[iprof].cs), sizeof(std::complex<double>));
    float HSBotrho = static_cast<float>(params.HSBot[iprof].rho); // 转换double为float
    MODFile.write(reinterpret_cast<char *>(&HSBotrho), sizeof(float));
    float SSPzn = static_cast<float>(params.SSP[iprof].z(params.SSP[iprof].z.size() - 1)); // 转换double为float
    MODFile.write(reinterpret_cast<char *>(&SSPzn), sizeof(float));
    for (int mode = 0; mode < eigen.M; mode++)
    {
        MODFile.seekp((iRecProf + 2 + mode) * 4 * LRecl, std::ios::beg);
        for (int irz = 0; irz < params.Pos->NRz; irz++)
        {
            std::complex<float> phiR = static_cast<std::complex<float>>(eigen.PsiR(mode, irz)); // 转换double为float
            MODFile.write(reinterpret_cast<char *>(&phiR), sizeof(std::complex<float>));
        }
    }

    // 写入模态数 M
    MODFile.seekp(iRecProf * 4 * LRecl, std::ios::beg);
    MODFile.write(reinterpret_cast<const char *>(&eigen.M), sizeof(int));

    // 写入复本征值 k
    int IFirst = 0; // C++ 从 0 开始
    for (int IREC = 0; IREC < (2 * eigen.M - 1) / LRecl + 1; ++IREC)
    {
        int ILast = std::min(int(eigen.M), IFirst + LRecl / 2) - 1;
        int segLen = ILast - IFirst + 1;

        // 定位到对应记录
        MODFile.seekp((iRecProf + 2 + eigen.M + IREC) * 4 * LRecl, std::ios::beg);

        // 将 eigen.k 中 IFirst 到 ILast 的复数写出
        // 转换为complex float输出

        for (int i = 0; i < segLen; i++)
        {
            complex<float> kf;
            kf = std::complex<float>(eigen.k(IFirst + i).real(), eigen.k(IFirst + i).imag());
            MODFile.write(reinterpret_cast<const char *>(&kf),
                          sizeof(std::complex<float>));
        }

        IFirst = ILast + 1;
    }

    // 更新下一段起始记录号
    iRecProf += 3 + eigen.M + (2 * eigen.M - 1) / LRecl;

    MODFile.close();
}

void kkc_interface::export_shd(std::string filename, int dataType)
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
    std::vector<double> freqvec(Nfreq, input.freqinfo->freq);
    std::vector<float> theta(Ntheta, 0.0), Sx(NSx, 0.0), Sy(NSy, 0.0);
    std::vector<float> Sz(input.Pos->NSz);
    for (int i = 0; i < input.Pos->NSz; i++)
    {
        Sz[i] = input.Pos->Sz[i];
    }
    std::vector<float> Rz(input.Pos->NRz);
    for (int i = 0; i < input.Pos->NRz; i++)
    {
        Rz[i] = input.Pos->Rz[i];
    }
    std::vector<float> Rr(input.Pos->NRr);
    for (int i = 0; i < input.Pos->NRr; i++)
    {
        Rr[i] = input.Pos->Rr[i];
    }

    std::string plottype;
    switch (input.Pos->GridType)
    {
    case Grid_Mode::MODE_I_Irregular:
        plottype = "irregular ";
        break;
    case Grid_Mode::MODE_R_Rectangular:
        plottype = "rectilin  ";
        break;
    }

    LRecl = std::max(std::max(std::max(std::max(std::max(41, 2 * input.freqinfo->Nfreq), Ntheta), input.Pos->NSz), input.Pos->NRz), 2 * input.Pos->NRr);

    SHDFile.write(reinterpret_cast<char *>(&LRecl), sizeof(int));
    SHDFile.write(input.Title.c_str(), input.Title.size());
    SHDFile.seekp(1 * 4 * LRecl, std::ios::beg);
    SHDFile.write(plottype.c_str(), 80);
    SHDFile.seekp(2 * 4 * LRecl, std::ios::beg);
    SHDFile.write(reinterpret_cast<char *>(&Nfreq), sizeof(int));
    SHDFile.write(reinterpret_cast<char *>(&Ntheta), sizeof(int));
    SHDFile.write(reinterpret_cast<char *>(&NSx), sizeof(int));
    SHDFile.write(reinterpret_cast<char *>(&NSy), sizeof(int));
    SHDFile.write(reinterpret_cast<char *>(&input.Pos->NSz), sizeof(int));
    SHDFile.write(reinterpret_cast<char *>(&input.Pos->NRz), sizeof(int));
    SHDFile.write(reinterpret_cast<char *>(&input.Pos->NRr), sizeof(int));
    float freq = static_cast<float>(input.freqinfo->freq); // 转换double为float
    SHDFile.write(reinterpret_cast<char *>(&freq), sizeof(float));
    SHDFile.write(reinterpret_cast<char *>(&Atten), sizeof(float));
    SHDFile.seekp(3 * 4 * LRecl, std::ios::beg);
    SHDFile.write(reinterpret_cast<char *>(freqvec.data()), Nfreq * sizeof(double));
    SHDFile.seekp(4 * 4 * LRecl, std::ios::beg);
    SHDFile.write(reinterpret_cast<char *>(theta.data()), Ntheta * sizeof(float));
    SHDFile.seekp(5 * 4 * LRecl, std::ios::beg);
    SHDFile.write(reinterpret_cast<char *>(Sx.data()), input.Pos->NSx * sizeof(float));
    SHDFile.seekp(6 * 4 * LRecl, std::ios::beg);
    SHDFile.write(reinterpret_cast<char *>(Sy.data()), input.Pos->NSy * sizeof(float));
    SHDFile.seekp(7 * 4 * LRecl, std::ios::beg);
    SHDFile.write(reinterpret_cast<char *>(Sz.data()), input.Pos->NSz * sizeof(float));
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
    for (int isz = 0; isz < input.Pos->NSz; ++isz)
    {
        for (int irz = 0; irz < input.Pos->NRz_per_range; ++irz)
        {
            int recnum = 10 + isz * input.Pos->NRz_per_range + irz;
            SHDFile.seekp(recnum * 4 * LRecl, std::ios::beg);
            for (int ir = 0; ir < input.Pos->NRr; ++ir)
            {
                std::complex<float> &P = data[GetFieldAddr(isz, irz, ir, input.Pos.get())];
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

parameters &kkc_interface::getParams() // 获取参数的引用
{
    return *params;
} // 获取参数


const parameters &kkc_interface::getParams_const() const
{
    return *params;
} // 获取参数

kkc_output &kkc_interface::getOutput() const
{
    return *output;
} // 获取输出
kkc_output kkc_interface::getOutput_Copy() const // 获取输出的副本
{
    return *output;
}

const kkc_output &kkc_interface::getOutput_const() const
{
    return *output;
} // 获取输出

bool kkc_interface::from_json(const std::string &jsonPath) // 从json读取参数
{
    // 按照setXX来定义json
    bool is_Success = false;
    std::ifstream file(jsonPath); // 打开文件流
    kkc_json j;                   // 创建json结构
    if (!file.is_open())
    {
        throw std::runtime_error("无法打开文件: " + jsonPath);
    }
    kkc_json in;
    file >> in; // 读取 JSON 数据
    file.close();

    // 解析频率点
    if (in.contains("freq"))
    {
        this->set_Freq(in["freq"].get<double>());
    }
    // 解析标题
    if (in.contains("Title"))
    {
        std::string title = in["Title"].get<std::string>();
        this->set_Title(title);
    }

    // 解析是否计算振速
    if (in.contains("is_Velocity"))
    {
        this->set_Velocity_enable(in["is_Velocity"].get<bool>());
    }

    // 调用各个模块的逆向解析函数
    // json_to_positions(in);
    // json_to_ssp(in);
    // json_to_boundary(in);
    // json_to_refl(in);
    // json_to_SBP(in);

    is_Success = true;
    return is_Success;
}

bool kkc_interface::to_json(const std::string &jsonPath) const // 将参数写入json
{
    bool is_Success = false;
    std::string jsonStr = this->to_json_string(); // 调用to_json_string()函数生成json字符串

    std::ofstream file(jsonPath); // 打开文件流
    if (file.is_open())
    {
        file << jsonStr; // 写入
        file.close();    // 关闭文件流
    }
    return true;
}

std::string kkc_interface::to_json_string() const // 将参数写入json字符串
{
    auto &params = this->getParams_const(); // 获取参数的常量引用
    kkc_json out;
    out["freq"] = params.freqinfo->freq;                                         // 频率
    out["Title"] = params.Title;                                                 // 标题
    params.is_Velocity ? out["is_Velocity"] = true : out["is_Velocity"] = false; // 是否计算振速

    // positions_to_json(out); // 位置
    // boundary_to_json(out);  // 边界条件
    // ssp_to_json(out);       // SSP
    // refl_to_json(out);      // 反射系数
    // SBP_to_json(out);       // 指向性
    return out.dump(4); // 写入json字符串
}