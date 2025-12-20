// 参数的基类
#ifndef PARAMSBASE_H
#define PARAMSBASE_H
#include "kkc_params.h"
#include "AttenMod.h"
#include "pchipMod.h"
#include "sspMod.h"
#include "util.h"

class paramsBase
{
public:
    paramsBase() {}
    virtual ~paramsBase() {}

    // 初始化逻辑
    virtual void Init(parameters &) const {}

    // 预设置逻辑
    virtual void SetupPre(parameters &) const {}

    // 设置默认参数值  必须重写的纯虚函数
    virtual void Default(parameters &) const = 0;

    // 读取参数逻辑
    virtual void Read(parameters &) const {}

    // 写入参数逻辑
    virtual void Write(parameters &) const {}

    // 后设置逻辑
    virtual void SetupPost(parameters &) const {}

    // 参数验证逻辑
    virtual void Validate(parameters &) const {}

    // 回显参数信息
    virtual void Echo(parameters &) const {}

    // 参数预处理逻辑
    virtual void Preprocess(parameters &) const {}

    // 清理资源逻辑
    virtual void Finalize(parameters &) const {}
};

class outputBase
{
public:
    outputBase() {}
    virtual ~outputBase() {}

    virtual void Init(kkc_output &) const {}
    /// Preprocessing as part of run.
    virtual void Preprocess(parameters &, kkc_output &) const {}
    /// Run the simulation.
    virtual void Run(parameters &, kkc_output &) const {};
    /// Postprocess after run is complete.
    virtual void Postprocess(parameters &, kkc_output &) const {}
    // 清理结果但不删除内存
    virtual void ClearResults(parameters &, kkc_output &) const {}
    // 清理资源逻辑
    virtual void Finalize(kkc_output &) const {}
};

#endif
