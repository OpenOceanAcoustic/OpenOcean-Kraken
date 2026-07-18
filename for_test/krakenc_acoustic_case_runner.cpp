#include "algorithm/AcousticCase.h"
#include "algorithm/KrakencSolver.h"

#include <nlohmann/json.hpp>

#include <exception>
#include <filesystem>
#include <iostream>
#include <string>

int main(int argc, char **argv)
{
    if (argc != 2 && argc != 3)
    {
        std::cerr << "Usage: OpenOceanKrakenc_acoustic_case_runner <case.env> [--diagnostics]\n";
        return 2;
    }
    const bool includeDiagnostics = argc == 3 &&
        std::string(argv[2]) == "--diagnostics";
    if (argc == 3 && !includeDiagnostics)
    {
        std::cerr << "unknown acoustic runner option: " << argv[2] << '\n';
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
        if (includeDiagnostics)
        {
            output["mesh_mode_counts"] = result.meshModeCounts;
            output["mesh_wavenumber_sets"] = nlohmann::ordered_json::array();
            for (const auto &set : result.meshWavenumberSets)
            {
                nlohmann::ordered_json values = nlohmann::ordered_json::array();
                for (const std::complex<double> value : set)
                {
                    values.push_back({value.real(), value.imag()});
                }
                output["mesh_wavenumber_sets"].push_back(std::move(values));
            }
        }
        output["wavenumbers"] = nlohmann::ordered_json::array();
        if (includeDiagnostics)
        {
            output["mode_diagnostics"] = nlohmann::ordered_json::array();
        }
        for (const OpenOceanKrakenc::ModeRoot &mode : result.modes)
        {
            output["wavenumbers"].push_back(
                {mode.wavenumber.real(), mode.wavenumber.imag()});
            if (includeDiagnostics)
            {
                output["mode_diagnostics"].push_back({
                    {"mesh_multiplier", mode.meshMultiplier},
                    {"matched_mesh_sets", mode.matchedMeshSets},
                    {"search_interval", {mode.searchLowerWavenumber,
                                          mode.searchUpperWavenumber}},
                    {"log10_residual", mode.diagnostic.log10Residual},
                    {"relative_correction", mode.diagnostic.relativeCorrection},
                    {"mesh_relative_spread", mode.meshRelativeSpread},
                    {"decision", mode.acceptanceDiagnostic}});
            }
        }
        if (includeDiagnostics)
        {
            output["root_search_diagnostics"] = nlohmann::ordered_json::array();
            for (const OpenOceanKrakenc::RootSearchDiagnostic &diagnostic :
                 result.rootSearchDiagnostics)
            {
                output["root_search_diagnostics"].push_back({
                    {"mesh_multiplier", diagnostic.meshMultiplier},
                    {"initial_guess", {diagnostic.initialGuess.real(),
                                        diagnostic.initialGuess.imag()}},
                    {"candidate", {diagnostic.candidate.real(),
                                    diagnostic.candidate.imag()}},
                    {"log10_residual", diagnostic.log10Residual},
                    {"relative_correction", diagnostic.relativeCorrection},
                    {"search_interval", {diagnostic.searchLowerWavenumber,
                                          diagnostic.searchUpperWavenumber}},
                    {"accepted", diagnostic.accepted},
                    {"decision", diagnostic.decision}});
            }
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
