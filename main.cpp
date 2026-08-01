#include "algorithm/AcousticCase.h"
#include "algorithm/KrakencSolver.h"
#include "algorithm/ModeFileWriter.h"
#include "algorithm/FieldSolver.h"
#include "OpenOceanKrakencKernelInterface.h"
#include "module/ParameterAdapters.h"
#include "module/json_in_out.hpp"

#include <nlohmann/json.hpp>

#include <algorithm>
#include <cctype>
#include <exception>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{
const char *loadErrorCodeName(
    OpenOceanKrakenc::LoadErrorCode code) noexcept
{
    using OpenOceanKrakenc::LoadErrorCode;
    switch (code)
    {
    case LoadErrorCode::MissingFile:
        return "missing_file";
    case LoadErrorCode::IoError:
        return "io_error";
    case LoadErrorCode::ParseError:
        return "parse_error";
    case LoadErrorCode::SchemaError:
        return "schema_error";
    case LoadErrorCode::InvalidField:
        return "invalid_field";
    case LoadErrorCode::UnsupportedCapability:
        return "unsupported_capability";
    case LoadErrorCode::None:
        return "none";
    }
    return "unknown";
}

std::string formatLoadError(
    const OpenOceanKrakenc::LoadError &error)
{
    return std::string("[") +
           loadErrorCodeName(error.code) +
           "] path='" + error.path + "': " + error.message;
}

std::vector<OpenOceanKrakenc::AcousticCase> loadAcousticCases(
    const std::filesystem::path &path)
{
    if (path.extension() == ".env")
    {
        return OpenOceanKrakenc::readAcousticEnvironments(path);
    }
    if (path.extension() == ".json")
    {
        OpenOceanKrakenc::OOKC_parameters params;
        const OpenOceanKrakenc::LoadResult loaded =
            OpenOceanKrakenc::read_json_file_result(
                path.string(), params);
        if (!loaded.ok)
        {
            throw std::runtime_error(
                "unable to read Krakenc JSON input " +
                formatLoadError(loaded.error));
        }
        return OpenOceanKrakenc::toAcousticCases(params);
    }
    throw std::invalid_argument("input extension must be .env or .json");
}

bool acousticInputExtension(const std::filesystem::path &path)
{
    return path.extension() == ".env" || path.extension() == ".json";
}
}

int main(int argc, char **argv)
{
    if (argc == 2 && std::string(argv[1]) == "--version")
    {
        std::cout << "OpenOcean-Krakenc 0.1.0\n";
        return 0;
    }
    if (argc == 1 || (argc == 2 && std::string(argv[1]) == "--help"))
    {
        std::cout << "Usage:\n"
                     "  OpenOcean-Krakenc --version\n"
                     "  OpenOcean-Krakenc --eigen <case.env|case.json>\n"
                     "  OpenOcean-Krakenc --mod <case.env|case.json> <output.mod> [--threads N]\n"
                     "  OpenOcean-Krakenc --field <case-root|case.env|case.json>\n"
                     "  OpenOcean-Krakenc --convert-env <case.env> <output.json>\n";
        return 0;
    }
    if (argc == 4 && std::string(argv[1]) == "--convert-env")
    {
        if (std::filesystem::path(argv[2]).extension() != ".env" ||
            std::filesystem::path(argv[3]).extension() != ".json")
        {
            std::cerr << "OpenOcean-Krakenc: --convert-env requires .env and .json paths\n";
            return 2;
        }
        try
        {
            OpenOceanKrakenc::KernelInterface api;
            const OpenOceanKrakenc::LoadResult loaded = api.loadEnv(argv[2]);
            if (!loaded.ok)
            {
                throw std::runtime_error("ENV conversion failed: " +
                                         loaded.error.message);
            }
            if (!api.to_json(argv[3]))
            {
                throw std::runtime_error("ENV conversion output could not be written");
            }
            nlohmann::ordered_json output;
            output["passes"] = true;
            output["path"] = std::filesystem::absolute(argv[3]).string();
            output["profile_count"] = api.getParams_const().sspInput.size();
            std::cout << output.dump() << '\n';
            return 0;
        }
        catch (const std::exception &error)
        {
            std::cerr << "OpenOcean-Krakenc ENV conversion failed: " << error.what() << '\n';
            return 1;
        }
    }
    if ((argc == 4 || argc == 6) && std::string(argv[1]) == "--mod")
    {
        if (!acousticInputExtension(argv[2]))
        {
            std::cerr << "OpenOcean-Krakenc: MOD input extension must be .env or .json\n";
            return 2;
        }
        try
        {
            std::size_t threadCount = 1;
            if (argc == 6)
            {
                if (std::string(argv[4]) != "--threads")
                {
                    throw std::invalid_argument("expected --threads N after MOD output");
                }
                const std::string value = argv[5];
                if (value.empty() || !std::all_of(
                        value.begin(), value.end(), [](unsigned char character) {
                            return std::isdigit(character) != 0;
                        }))
                {
                    throw std::invalid_argument("MOD thread count must be a positive integer");
                }
                threadCount = std::stoull(value);
                if (threadCount == 0 || threadCount > 1024)
                {
                    throw std::invalid_argument("MOD thread count must be between 1 and 1024");
                }
            }
            const std::vector<OpenOceanKrakenc::AcousticCase> inputs =
                loadAcousticCases(std::filesystem::path(argv[2]));
            const std::vector<OpenOceanKrakenc::ModeFileWriteResult> results =
                OpenOceanKrakenc::writeModeFile(
                    inputs, std::filesystem::path(argv[3]), threadCount);
            nlohmann::ordered_json output;
            output["passes"] = true;
            output["profile_count"] = results.size();
            output["mode_count"] = results.front().modeCount;
            output["tabulated_depth_count"] = results.front().tabulatedDepthCount;
            output["record_length_words"] = results.front().recordLengthWords;
            output["worker_threads"] = std::min(threadCount, inputs.size());
            output["requested_threads"] = threadCount;
            output["profiles"] = nlohmann::ordered_json::array();
            for (const auto &result : results)
            {
                output["profiles"].push_back({
                    {"mode_count", result.modeCount},
                    {"tabulated_depth_count", result.tabulatedDepthCount},
                    {"record_length_words", result.recordLengthWords}});
            }
            output["path"] = std::filesystem::absolute(argv[3]).string();
            std::cout << output.dump() << '\n';
            return 0;
        }
        catch (const std::exception &error)
        {
            std::cerr << "OpenOcean-Krakenc MOD write failed: " << error.what() << '\n';
            return 1;
        }
    }
    if (argc == 3 && std::string(argv[1]) == "--eigen")
    {
        if (!acousticInputExtension(argv[2]))
        {
            std::cerr << "OpenOcean-Krakenc: eigen input extension must be .env or .json\n";
            return 2;
        }
        try
        {
            const std::vector<OpenOceanKrakenc::AcousticCase> inputs =
                loadAcousticCases(std::filesystem::path(argv[2]));
            const OpenOceanKrakenc::AcousticCase &input = inputs.front();
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
            std::cerr << "OpenOcean-Krakenc eigen solve failed: " << error.what() << '\n';
            return 1;
        }
    }
    if (argc == 3 && std::string(argv[1]) == "--field")
    {
        try
        {
            const std::filesystem::path inputPath = argv[2];
            if (inputPath.extension() == ".env" || inputPath.extension() == ".json")
            {
                OpenOceanKrakenc::KernelInterface api;
                const OpenOceanKrakenc::LoadResult loaded =
                    inputPath.extension() == ".env"
                        ? api.loadEnv(inputPath.string())
                        : api.loadJson(inputPath.string());
                if (!loaded.ok)
                {
                    throw std::runtime_error(
                        "unable to read Krakenc ENV/JSON field input: " +
                        loaded.error.message);
                }
                api.runField();
                std::filesystem::path outputPath = inputPath;
                outputPath.replace_extension(".shd");
                api.export_shd(outputPath.string(), 1);
                nlohmann::ordered_json output;
                output["passes"] = true;
                output["source_depth_count"] = api.getOutput_const().sourceCount;
                output["receiver_depth_count"] = api.getOutput_const().receiverDepthCount;
                output["range_count"] = api.getOutput_const().rangeCount;
                output["path"] = std::filesystem::absolute(outputPath).string();
                std::cout << output.dump() << '\n';
                return 0;
            }
            if (!inputPath.extension().empty())
            {
                std::cerr << "OpenOcean-Krakenc: field input must be a root path, .env, or .json\n";
                return 2;
            }
            const std::string root = argv[2];
            const OpenOceanKrakenc::ModeFileData modes =
                OpenOceanKrakenc::readModeFile(root + ".mod");
            OpenOceanKrakenc::FieldParameters parameters =
                OpenOceanKrakenc::readFieldParameters(root + ".flp");
            if (std::filesystem::exists(root + ".env"))
            {
                const OpenOceanKrakenc::AcousticCase sourceProfile =
                    OpenOceanKrakenc::readAcousticEnv(root + ".env");
                for (double depth : parameters.sourceDepths)
                {
                    parameters.sourceSoundSpeeds.push_back(
                        OpenOceanKrakenc::soundSpeedAt(sourceProfile, depth));
                }
            }
            const OpenOceanKrakenc::PressureField field =
                OpenOceanKrakenc::evaluateField(modes, parameters);
            OpenOceanKrakenc::writeShadeFile(field, root + ".shd");
            nlohmann::ordered_json output;
            output["passes"] = true;
            output["source_depth_count"] = field.sourceDepths.size();
            output["receiver_depth_count"] = field.receiverDepths.size();
            output["range_count"] = field.rangesMetres.size();
            output["path"] = std::filesystem::absolute(root + ".shd").string();
            std::cout << output.dump() << '\n';
            return 0;
        }
        catch (const std::exception &error)
        {
            std::cerr << "OpenOcean-Krakenc field solve failed: " << error.what() << '\n';
            return 1;
        }
    }
    std::cerr << "OpenOcean-Krakenc: unsupported arguments; use --help\n";
    return 2;
}
