#include "algorithm/AcousticCase.h"
#include "algorithm/ComplexDispersion.h"
#include "algorithm/ComplexMatrixBuilder.h"
#include "algorithm/ComplexRootFinder.h"
#include "algorithm/FieldSolver.h"

#include <nlohmann/json.hpp>

#include <complex>
#include <filesystem>
#include <iostream>
#include <vector>

int main(int argc, char **argv)
{
    if (argc != 3)
    {
        std::cerr << "Usage: OpenOceanKrakenc_root_audit_runner <case.env> <oracle.mod>\n";
        return 2;
    }
    try
    {
        using namespace OpenOceanKrakenc;
        const AcousticCase input = readAcousticEnv(argv[1]);
        const ModeFileData oracle = readModeFile(argv[2]);
        nlohmann::ordered_json output;
        output["passes"] = true;
        output["mesh_audits"] = nlohmann::ordered_json::array();
        for (const int multiplier : {1, 2})
        {
            const AcousticMatrix matrix = buildAcousticMatrix(input, multiplier);
            const ComplexDispersion dispersion(input, matrix);
            std::vector<std::complex<double>> accepted;
            nlohmann::ordered_json roots = nlohmann::ordered_json::array();
            for (const std::complex<double> expectedK : oracle.wavenumbers)
            {
                const std::complex<double> guess = expectedK * expectedK;
                const double tolerance = std::abs(guess) *
                    static_cast<double>(matrix.b1.size()) * 1.0e-14;
                const RootResult root = complexSecantDispersion(
                    guess, tolerance, 1000, dispersion, accepted);
                const std::complex<double> actualK = std::sqrt(root.root);
                roots.push_back({
                    {"expected", {expectedK.real(), expectedK.imag()}},
                    {"actual", {actualK.real(), actualK.imag()}},
                    {"converged", root.converged},
                    {"log10_residual", root.log10Residual},
                    {"relative_difference", std::abs(actualK - expectedK) /
                        std::max(1.0e-30, std::abs(expectedK))}});
                if (root.converged)
                {
                    accepted.push_back(root.root);
                }
            }
            output["mesh_audits"].push_back({
                {"multiplier", multiplier}, {"roots", std::move(roots)}});
        }
        std::cout << output.dump() << '\n';
        return 0;
    }
    catch (const std::exception &error)
    {
        std::cerr << "root audit failed: " << error.what() << '\n';
        return 1;
    }
}
