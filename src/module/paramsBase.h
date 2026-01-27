// 参数的基类
#ifndef PARAMSBASE_H
#define PARAMSBASE_H
#include "OpenOceanKrakenParams.h"
#include "AttenMod.h"
#include "pchipMod.h"
#include "sspMod.h"
#include "util.h"

namespace OpenOceanKraken
{
    class paramsBase
    {
    public:
        paramsBase() {}
        virtual ~paramsBase() {}

        // 初始化逻辑
        virtual void Init(OOK_parameters &) const {}

        // 预设置逻辑
        virtual void SetupPre(OOK_parameters &) const {}

        // 设置默认参数值  必须重写的纯虚函数
        virtual void Default(OOK_parameters &) const = 0;

        // 读取参数逻辑
        virtual void Read(OOK_parameters &) const {}

        // 写入参数逻辑
        virtual void Write(OOK_parameters &) const {}

        // 后设置逻辑
        virtual void SetupPost(OOK_parameters &) const {}

        // 参数验证逻辑
        virtual void Validate(OOK_parameters &) const {}

        // 回显参数信息
        virtual void Echo(OOK_parameters &) const {}

        // 参数预处理逻辑
        virtual void Preprocess(OOK_parameters &) const {}

        // 清理资源逻辑
        virtual void Finalize(OOK_parameters &) const {}
    };

    class outputBase
    {
    public:
        outputBase() {}
        virtual ~outputBase() {}

        virtual void Init(OOK_output &) const {}
        /// Preprocessing as part of run.
        virtual void Preprocess(OOK_parameters &, OOK_output &) const {}
        /// Run the simulation.
        virtual void Run(OOK_parameters &, OOK_output &) const {};
        /// Postprocess after run is complete.
        virtual void Postprocess(OOK_parameters &, OOK_output &) const {}
        // 清理结果但不删除内存
        virtual void ClearResults(OOK_parameters &, OOK_output &) const {}
        // 清理资源逻辑
        virtual void Finalize(OOK_output &) const {}
    };

}
#endif
