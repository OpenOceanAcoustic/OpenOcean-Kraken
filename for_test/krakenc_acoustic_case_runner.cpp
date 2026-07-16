#include "algorithm/AcousticCase.h"
#include "algorithm/KrakencSolver.h"

#include <nlohmann/json.hpp>

#include <exception>
#include <filesystem>
#include <iostream>

int main(int argc, char **argv)
{
    if (argc != 2)
    {
        std::cerr << "Usage: OpenOceanKrakenc_acoustic_case_runner <case.env>\n";
        return 2;
    }

    try
    {
        const OpenOceanKrakenc::AcousticCase input =
            OpenOceanKrakenc::readAcousticEnv(std::filesystem::path(argv[1]));
        const OpenOceanKrakenc::AcousticSolveResult result =
            OpenOceanKrakenc::solveAcousticModes(input);

        nlohmann::ordered_json output;
        output["passes"] = true;
        output["mode_count"] = result.modes.size();
        output["mesh_sets_used"] = result.meshSetsUsed;
        output["wavenumbers"] = nlohmann::ordered_json::array();
        for (const OpenOceanKrakenc::ModeRoot &mode : result.modes)
        {
            output["wavenumbers"].push_back(
                {mode.wavenumber.real(), mode.wavenumber.imag()});
        }
        std::cout << output.dump() << '\n';
        return 0;
    }
    catch (const std::exception &error)
    {
        std::cerr << "OpenOcean-Krakenc acoustic solve failed: " << error.what() << '\n';
        return 1;
    }
}
