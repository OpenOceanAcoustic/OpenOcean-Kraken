#include "OpenOceanKrakencInterface.h"
#include "algorithm/FieldSolver.h"

#include <cmath>
#include <complex>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace OpenOceanKrakenc;

namespace
{
void require(bool condition, const char *message)
{
    if (!condition) throw std::runtime_error(message);
}

template <typename Value>
Value readValue(const std::vector<char> &data, std::size_t offset)
{
    Value value{};
    std::memcpy(&value, data.data() + offset, sizeof(Value));
    return value;
}

std::vector<std::complex<float>> readShade(const std::filesystem::path &path,
                                           std::string &plotType)
{
    std::ifstream stream(path, std::ios::binary);
    require(static_cast<bool>(stream), "unable to open exported SHD");
    std::vector<char> data((std::istreambuf_iterator<char>(stream)), {});
    const int words = readValue<int>(data, 0);
    const std::size_t bytes = static_cast<std::size_t>(4 * words);
    plotType.assign(data.data() + bytes, data.data() + bytes + 10);
    const std::size_t header = 2 * bytes;
    const int nsz = readValue<int>(data, header + 16);
    const int nrz = readValue<int>(data, header + 20);
    const int nrr = readValue<int>(data, header + 24);
    std::vector<std::complex<float>> values;
    values.reserve(static_cast<std::size_t>(nsz * nrz * nrr));
    std::size_t record = 10;
    for (int source = 0; source < nsz; ++source)
    {
        for (int receiver = 0; receiver < nrz; ++receiver, ++record)
        {
            for (int range = 0; range < nrr; ++range)
            {
                const std::size_t offset = record * bytes + static_cast<std::size_t>(8 * range);
                values.emplace_back(readValue<float>(data, offset),
                                    readValue<float>(data, offset + 4));
            }
        }
    }
    return values;
}
}

int main()
{
    const std::filesystem::path source = OPENOCEANKRAKENC_SOURCE_DIR;
    Interface api;
    require(api.from_env((source / "for_test/fixtures/two_profile_small.env").string()),
            "export fixture ENV load failed");
    api.getParams().shdPath.clear();
    api.set_Velocity_enable(true);
    api.run();
    api.getParams().envPath.clear();
    api.getParams().flpPath.clear();

    const std::filesystem::path directory =
        std::filesystem::temp_directory_path() / "openocean_krakenc_export_test";
    std::filesystem::create_directories(directory);
    const std::filesystem::path modeRoot = directory / "memory_mode";
    api.export_mod(modeRoot.string());
    const ModeFileData modes = readModeFile(modeRoot.string() + ".mod");
    require(modes.wavenumbers.size() == static_cast<std::size_t>(api.getOutput_const().eigen.front().M),
            "exported MOD mode count mismatch");
    require(modes.additionalProfiles.size() == 1, "exported MOD profile count mismatch");

    const char *names[] = {"pressure", "vertical", "horizontal"};
    const std::vector<std::complex<float>> *expected[] = {
        &api.getOutput_const().pressure,
        &api.getOutput_const().verticalVelocity,
        &api.getOutput_const().horizontalVelocity};
    for (int type = 1; type <= 3; ++type)
    {
        const std::filesystem::path root = directory / names[type - 1];
        api.export_shd(root.string(), type);
        std::string plotType;
        const auto actual = readShade(root.string() + ".shd", plotType);
        require(plotType.find("rectilin") != std::string::npos,
                "SHD plot type was not written");
        require(actual == *expected[type - 1], "SHD data does not match memory output");
    }

    bool rejected = false;
    try { api.export_shd((directory / "bad").string(), 4); }
    catch (const std::invalid_argument &) { rejected = true; }
    require(rejected, "invalid SHD data type was accepted");

    const std::filesystem::path combined = directory / "combined";
    api.export_result(combined.string());
    require(std::filesystem::exists(combined.string() + "_P.mod"), "combined MOD missing");
    require(std::filesystem::exists(combined.string() + "_P.shd"), "combined pressure SHD missing");
    require(std::filesystem::exists(combined.string() + "_V.shd"), "combined vertical SHD missing");
    require(std::filesystem::exists(combined.string() + "_H.shd"), "combined horizontal SHD missing");

    const std::filesystem::path dotted = directory / "combined.v1";
    api.export_result(dotted.string());
    require(std::filesystem::exists(dotted.string() + "_P.mod"),
            "dotted combined MOD missing");
    require(std::filesystem::exists(dotted.string() + "_P.shd"),
            "dotted combined pressure SHD missing");
    require(std::filesystem::exists(dotted.string() + "_V.shd"),
            "dotted combined vertical SHD missing");
    require(std::filesystem::exists(dotted.string() + "_H.shd"),
            "dotted combined horizontal SHD missing");

    Interface pressureOnly;
    require(pressureOnly.from_env((source / "for_test/fixtures/two_profile_small.env").string()),
            "pressure-only export ENV load failed");
    pressureOnly.getParams().shdPath.clear();
    pressureOnly.runField();
    pressureOnly.export_shd((directory / "pressure_only").string(), 1);
    bool missingVelocityRejected = false;
    try { pressureOnly.export_shd((directory / "missing_velocity").string(), 2); }
    catch (const std::logic_error &) { missingVelocityRejected = true; }
    require(missingVelocityRejected,
            "velocity SHD export succeeded while velocity calculation was disabled");

    Interface elastic;
    const std::filesystem::path elasticEnv =
        source.parent_path() / "test" / "elastic_fd_two_layer.env";
    require(elastic.from_env(elasticEnv.string()), "elastic export ENV load failed");
    elastic.runEigen();
    const std::filesystem::path elasticRoot = directory / "elastic_mode";
    elastic.export_mod(elasticRoot.string());
    const ModeFileData elasticModes =
        readModeFile(elasticRoot.string() + ".mod");
    require(elasticModes.meshCounts.size() == 1,
            "elastic MOD export must describe only the acoustic mode media");
    require(elasticModes.materials.size() == 1 &&
                elasticModes.materials.front() == "ACOUSTIC",
            "elastic MOD export wrote non-acoustic mode media");
    require(std::abs(elasticModes.top.depth) < 1.0e-6 &&
                std::abs(elasticModes.bottom.depth - 200.0) < 1.0e-6,
            "elastic MOD halfspace depths do not bracket the complete media stack");

    std::error_code ignored;
    std::filesystem::remove_all(directory, ignored);
    return 0;
}
