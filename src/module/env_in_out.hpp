#ifndef ENV_IN_OUT_HPP
#define ENV_IN_OUT_HPP

#include "OpenOceanKrakenParams.h"

namespace OpenOceanKraken
{
    // 读取 env 文件并填充 OOK_parameters
    bool read_env_file(const std::string &envPath, OOK_parameters &params);

}

#endif
