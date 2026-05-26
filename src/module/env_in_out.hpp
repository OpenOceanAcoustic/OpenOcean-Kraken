#ifndef ENV_IN_OUT_HPP
#define ENV_IN_OUT_HPP

#include "OpenOceanKrakenParams.h"

namespace OpenOceanKraken
{
    // 读取 sbp 文件（声源指向性）
    bool read_sbp_file(const std::string &envPath, OOK_parameters &params);

    // 读取 flp 文件（声源/接收器/运行选项）
    bool read_flp_file(const std::string &envPath, OOK_parameters &params);

    // 读取反射系数文件（.trc / .brc）
    bool read_refCoef_file(const std::string &envPath, OOK_parameters &params, std::string pattern);

    // 读取 env 文件并填充 OOK_parameters
    bool read_env_file(const std::string &envPath, OOK_parameters &params);

//     // 从 env 文件创建 OOK_parameters 对象（便捷函数）
//     OOK_parameters env_to_params(const std::string &envPath);
}

#endif
