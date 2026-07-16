#include "module/env_in_out.hpp"

#include "algorithm/AcousticCase.h"
#include "algorithm/FieldSolver.h"
#include "module/ParameterAdapters.h"

#include <filesystem>
#include <utility>

namespace OpenOceanKrakenc
{
bool read_env_file(const std::string &envPath, OOKC_parameters &params)
{
    try
    {
        OOKC_parameters candidate = params;
        const std::vector<AcousticCase> cases = readAcousticEnvironments(envPath);
        updatePublicParameters(cases, candidate);
        candidate.envPath = std::filesystem::absolute(envPath).lexically_normal().string();
        std::filesystem::path root = candidate.envPath;
        root.replace_extension();
        candidate.flpPath = root.string() + ".flp";
        candidate.modPath = root.string() + ".mod";
        candidate.shdPath = root.string() + ".shd";
        params = std::move(candidate);
        return true;
    }
    catch (...)
    {
        return false;
    }
}

bool read_flp_file(const std::string &envPath, OOKC_parameters &params)
{
    try
    {
        std::filesystem::path flpPath = envPath;
        if (flpPath.extension() != ".flp")
        {
            flpPath.replace_extension(".flp");
        }
        OOKC_parameters candidate = params;
        const FieldParameters field = readFieldParameters(flpPath);
        updatePublicParameters(field, candidate);
        candidate.flpPath = std::filesystem::absolute(flpPath).lexically_normal().string();
        params = std::move(candidate);
        return true;
    }
    catch (...)
    {
        return false;
    }
}
}
