#ifndef ENV_IN_OUT_HPP
#define ENV_IN_OUT_HPP

#include "OpenOceanKrakenParams.h"

namespace OpenOceanKraken
{
    // 解析env文件并填充OOK_parameters结构
    bool env_to_json(const std::string &envPath, OOK_parameters &params);
    
    // 将env文件转换为JSON字符串
    std::string env_to_json_string(const std::string &envPath);
    
    // 将env文件转换为OpenOcean_json对象
    OpenOcean_json env_to_json_obj(const std::string &envPath);
    
    // 从env文件创建OOK_parameters对象
    OOK_parameters env_to_params(const std::string &envPath);
}

#endif