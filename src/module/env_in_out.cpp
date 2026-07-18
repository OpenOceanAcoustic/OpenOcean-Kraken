#include "module/env_in_out.hpp"

#include "algorithm/AcousticCase.h"
#include "algorithm/FieldSolver.h"
#include "module/ParameterAdapters.h"

#include <filesystem>
#include <utility>

namespace OpenOceanKrakenc
{
namespace
{
bool isRegularFile(const std::filesystem::path &path)
{
    std::error_code error;
    return std::filesystem::is_regular_file(path, error);
}

LoadResult inputFailure(const std::string &path, const std::exception &error)
{
    const std::string message = error.what();
    const LoadErrorCode code =
        message.find("invalid") != std::string::npos ||
                message.find("must") != std::string::npos
            ? LoadErrorCode::InvalidField
            : LoadErrorCode::ParseError;
    return LoadResult::failure(code, path, message);
}
}

LoadResult read_env_file_result(const std::string &envPath,
                                OOKC_parameters &params)
{
    const std::filesystem::path path = envPath;
    if (!isRegularFile(path))
    {
        return LoadResult::failure(LoadErrorCode::MissingFile, envPath,
                                   "ENV file does not exist or is not a regular file");
    }
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
        return LoadResult::success();
    }
    catch (const std::exception &error)
    {
        return inputFailure(envPath, error);
    }
}

LoadResult read_flp_file_result(const std::string &envPath,
                                OOKC_parameters &params)
{
    std::filesystem::path flpPath = envPath;
    if (flpPath.extension() != ".flp")
    {
        flpPath.replace_extension(".flp");
    }
    if (!isRegularFile(flpPath))
    {
        return LoadResult::failure(LoadErrorCode::MissingFile, flpPath.string(),
                                   "FLP file does not exist or is not a regular file");
    }
    try
    {
        OOKC_parameters candidate = params;
        const FieldParameters field = readFieldParameters(flpPath);
        updatePublicParameters(field, candidate);
        candidate.flpPath = std::filesystem::absolute(flpPath).lexically_normal().string();
        params = std::move(candidate);
        return LoadResult::success();
    }
    catch (const std::exception &error)
    {
        return inputFailure(flpPath.string(), error);
    }
}

bool read_env_file(const std::string &envPath, OOKC_parameters &params)
{
    return read_env_file_result(envPath, params).ok;
}

bool read_flp_file(const std::string &envPath, OOKC_parameters &params)
{
    return read_flp_file_result(envPath, params).ok;
}
}
