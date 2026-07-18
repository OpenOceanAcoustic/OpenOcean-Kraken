#include "OpenOceanKrakencInterface.h"
#include "module/env_in_out.hpp"

#include <cassert>
#include <filesystem>
#include <stdexcept>

using namespace OpenOceanKrakenc;

#undef assert
#define assert(condition) do { if (!(condition)) throw std::runtime_error("assertion failed: " #condition); } while (false)

int main()
{
    const std::filesystem::path workspace = OPENOCEANKRAKENC_WORKSPACE_DIR;
    const std::filesystem::path env =
        workspace / "test" / "toolbox_env" / "MunkKleaky.env";

    OOKC_parameters params;
    assert(read_env_file(env.string(), params));
    assert(std::filesystem::equivalent(params.envPath, env));
    assert(params.Title.find("Munk") != std::string::npos);
    assert(params.freqinfo.freq == 50.0);
    assert(params.NProf == 1);
    assert(params.sspInput.size() == 1);
    assert(params.Pos.NSz >= 1);
    assert(params.Pos.NRz >= 1);

    assert(read_flp_file(env.string(), params));
    assert(params.Pos.NRr == 501);
    assert(params.Pos.NSz == 2);
    assert(params.Pos.NRz == 1001);
    assert(params.SourceType == Source_Mode::MODE_R_Point);
    assert(params.modeType == ModeType::Adiabatic);
    assert(params.SBP.isSet);

    const OOKC_parameters before = params;
    assert(!read_env_file((workspace / "missing.env").string(), params));
    assert(params.Title == before.Title);
    assert(params.freqinfo.freq == before.freqinfo.freq);
    assert(params.sspInput.size() == before.sspInput.size());

    Interface api;
    assert(api.from_env(env.string()));
    assert(api.getParams_const().Pos.NRr == 501);

    const std::filesystem::path isolatedEnv =
        std::filesystem::temp_directory_path() / "openocean_krakenc_without_flp.env";
    std::filesystem::copy_file(env, isolatedEnv,
                               std::filesystem::copy_options::overwrite_existing);
    Interface missingFlp;
    missingFlp.set_Title("transaction sentinel");
    assert(!missingFlp.from_env(isolatedEnv.string()));
    assert(missingFlp.getParams_const().Title == "transaction sentinel");
    std::error_code ignored;
    std::filesystem::remove(isolatedEnv, ignored);
    return 0;
}
