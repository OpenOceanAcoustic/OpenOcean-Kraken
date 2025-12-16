#ifndef KKCIFACE_H
#define KKCIFACE_H

#include "kkc_params.h"
#include <memory>

class kkc_interface_PIMPL;

// 只作为一个参数表，不涉及具体的计算过程
class kkc_interface
{
public:
    kkc_interface();  // 构造函数
    ~kkc_interface(); // 析构函数
    void init(); // 初始化
    void setup(); // 配置
    void input_setup(); // 设置输入参数
    void output_setup(); // 设置输出参数
    void run(); // 运行
    void clearResults(); // 清除结果
    void runField(); // 运行声场
    void runSolveV(); // 运行特征值求解器
    void free(); // 释放内存
private:
    kkc_output *output; // 输出
    parameters *params; // 输入
    // std::atomic<int> sharedJobID;//共享任务ID

    int NumThreads;          // 线程数
    size_t totalTasks;       // 任务数
    size_t field_size;       // 场大小

    bool is_setup = false; // 是否设置了参数
    std::unique_ptr<kkc_interface_PIMPL> impl;  
};

#endif // KKCIFACE_H
