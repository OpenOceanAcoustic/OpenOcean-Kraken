#include "kkc_interface.h"

void print_header()
{
    std::cout << " - OpenOcean-Kraken: A parallel underwater acoustic nomal modes simulator with pybind11 Python and C++ interfaces\n"
                 "\n"
                 "Copyright (C) 2024-2025 OpenOceanAcoustic\n"
                 "Authors: Liang Yi, Qian Peng, Dai Nuoge, Luo Shixiong, Liu Ruihang, Liu Jiatong, Wang Zhengwei\n"
                 "Mail: yingxinliang1@gmail.com\n"
                 "Based on Kraken, which is Copyright (C) 1983-2024 Michael B. Porter\n"
                 "GPL3 licensed, no warranty, see LICENSE or https://www.gnu.org/licenses/\n"
                 "\n";
}

// 定义实现类
class kkc_interface_PIMPL
{
public:
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
}

// 配置
void kkc_interface::setup() // 设置参数
{
    input_setup(); // 设置输入参数(外界输入参数，因此每次计算时候都重新计算一次环境)
    if (!this->is_setup)
    {
        this->output_setup();  // 若未配置，则配置输出的内存
        this->is_setup = true; // 设置配置标志位
    }
}

// 输入配置
void kkc_interface::input_setup() // 设置输入参数
{
}

// 输出配置
void kkc_interface::output_setup() // 设置输出参数
{
}

void kkc_interface::run() // 运行
{
    print_header();             // 打印头信息
    this->setup();              // 配置
    this->runSolveV();      // 求解本征值和本征函数
    this->runField();           // 求解声场
    this->output_postprocess(); // 输出后处理
}

void kkc_interface::clearResults() // 清除结果
{
}

void kkc_interface::free() // 释放内存
{
}

void kkc_interface::runSolveV()
{
    EigenVWorker(...);
}

void kkc_interface::runField()
{
    FieldSolveWorker(...);
}