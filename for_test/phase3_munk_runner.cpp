#include "algorithm/AcousticCase.h"
#include "algorithm/ComplexMatrixBuilder.h"
#include "algorithm/ComplexModeNormalization.h"
#include "algorithm/ComplexModeSolver.h"
#include "algorithm/KrakencSolver.h"

#include <nlohmann/json.hpp>

#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>

using namespace OpenOceanKrakenc;

int main(int argc, char **argv)
{
    try
    {
        if (argc != 3)
        {
            throw std::invalid_argument("usage: phase3_munk_runner ENV OUTPUT_JSON");
        }
        const AcousticCase input = readAcousticEnv(argv[1]);
        const AcousticMatrix matrix = buildAcousticMatrix(input, 1);
        const AcousticSolveResult roots = solveAcousticModes(input);

        constexpr int sampleStride = 5;
        nlohmann::json output;
        output["mode_count"] = roots.modes.size();
        output["depths"] = nlohmann::json::array();
        for (std::size_t index = 0; index < matrix.depth.size(); index += sampleStride)
        {
            output["depths"].push_back(matrix.depth[index]);
        }
        output["modes"] = nlohmann::json::array();
        output["group_velocities"] = nlohmann::json::array();

        for (const ModeRoot &root : roots.modes)
        {
            const ComplexModeResult raw =
                solveAcousticMode(input, matrix, root.eigenvalue);
            if (!raw.converged)
            {
                throw std::runtime_error("complex mode inverse iteration failed");
            }
            const NormalizedComplexMode normalized = normalizeAcousticMode(
                input, matrix, root.eigenvalue, raw.turningPoint, raw.mode);
            nlohmann::json samples = nlohmann::json::array();
            for (Eigen::Index index = 0; index < normalized.mode.size();
                 index += sampleStride)
            {
                samples.push_back({normalized.mode[index].real(),
                                   normalized.mode[index].imag()});
            }
            output["modes"].push_back(std::move(samples));
            output["group_velocities"].push_back(normalized.groupVelocity);
        }

        std::ofstream stream(argv[2]);
        if (!stream)
        {
            throw std::runtime_error("cannot open phase-3 JSON output");
        }
        stream << output;
        return 0;
    }
    catch (const std::exception &error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
