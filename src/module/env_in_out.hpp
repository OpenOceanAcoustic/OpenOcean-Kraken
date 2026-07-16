#ifndef OPEN_OCEAN_KRAKENC_ENV_IN_OUT_HPP
#define OPEN_OCEAN_KRAKENC_ENV_IN_OUT_HPP

#include "OpenOceanKrakencParams.h"

#include <string>

namespace OpenOceanKrakenc
{
bool read_env_file(const std::string &envPath, OOKC_parameters &params);
bool read_flp_file(const std::string &envPath, OOKC_parameters &params);
}

#endif
