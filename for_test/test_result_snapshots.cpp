#include "OpenOceanKrakenKernelInterface.h"
#include "ThreadPool.h"

#include <algorithm>
#include <array>
#include <complex>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

using OpenOceanKraken::ThreadPool;

namespace
{
    void require(bool condition, const std::string &message)
    {
        if (!condition)
        {
            throw std::runtime_error(message);
        }
    }
}

int main(int argc, char **argv)
{
    if (argc != 2)
    {
        std::cerr << "Usage: OpenOceanKraken_snapshot_tests <fixture.env>\n";
        return 2;
    }

    try
    {
        ThreadPool pool(1);
        OpenOceanKraken::KernelInterface interface(pool);
        interface.setNumThreads(1);
        require(interface.from_env(argv[1]), "fixture must load");
        interface.set_Velocity_enable(true);
        interface.run();

        const auto pressure = interface.getPressureCopy();
        const auto &params = interface.getParams_const();
        require(pressure.source_count == static_cast<std::size_t>(params.Pos.NSz), "source count");
        require(pressure.range_count == static_cast<std::size_t>(params.Pos.NRr), "range count");
        require(pressure.depth_count == static_cast<std::size_t>(params.Pos.NRz_per_range), "depth count");
        require(
            pressure.values.size() == pressure.source_count * pressure.range_count * pressure.depth_count,
            "pressure value count");
        require(pressure.source_depths.size() == pressure.source_count, "source coordinates");
        require(pressure.receiver_ranges.size() == pressure.range_count, "range coordinates");
        require(pressure.receiver_depths.size() == pressure.depth_count, "depth coordinates");

        const auto vertical = interface.getVerticalVelocityCopy();
        const auto horizontal = interface.getHorizontalVelocityCopy();
        require(vertical.values.size() == pressure.values.size(), "vertical value count");
        require(horizontal.values.size() == pressure.values.size(), "horizontal value count");

        const auto modes = interface.getModesCopy();
        require(modes.size() == static_cast<std::size_t>(params.NProf), "profile count");
        require(!modes.front().wavenumbers.empty(), "mode wavenumbers");
        require(
            modes.front().mode_shapes.rows() == static_cast<Eigen::Index>(modes.front().wavenumbers.size()),
            "mode-shape row count");
        require(
            modes.front().mode_shapes.cols() == static_cast<Eigen::Index>(modes.front().depth.size()),
            "mode-shape depth count");

        const std::filesystem::path shadeRoot =
            std::filesystem::temp_directory_path() /
            "openocean_kraken_snapshot";
        interface.export_shd(shadeRoot.string(), 1);
        std::ifstream shade(
            shadeRoot.string() + ".shd",
            std::ios::binary);
        require(static_cast<bool>(shade), "SHD snapshot must open");
        int recordWords = 0;
        shade.read(
            reinterpret_cast<char *>(&recordWords),
            sizeof(recordWords));
        require(recordWords >= 20, "SHD record length");
        shade.seekg(
            static_cast<std::streamoff>(recordWords) * 4,
            std::ios::beg);
        std::array<char, 80> plotType{};
        shade.read(plotType.data(), plotType.size());
        require(
            std::string(plotType.data(), 10) == "rectilin  ",
            "SHD plot-type prefix");
        require(
            std::all_of(
                plotType.begin() + 10,
                plotType.end(),
                [](char value) { return value == '\0'; }),
            "SHD plot-type record must use deterministic zero padding");
        std::error_code ignored;
        std::filesystem::remove(shadeRoot.string() + ".shd", ignored);

        const std::complex<float> saved = pressure.values.front();
        interface.clearResults();
        require(pressure.values.front() == saved, "snapshot must own its values");

        interface.set_RProf(0.0, 10000.0, 3);
        interface.set_MLimit(64);
        interface.set_CoherenceType(OpenOceanKraken::CoherenceType::Coherent);
        interface.set_ModeType(OpenOceanKraken::ModeType::Adiabatic);
        require(interface.getParams_const().NProf == 3, "profile setter");
        require(interface.getParams_const().MLimit == 64, "mode limit setter");

        std::cout << "OOK snapshot tests passed\n";
        return 0;
    }
    catch (const std::exception &error)
    {
        std::cerr << "OOK snapshot test failed: " << error.what() << '\n';
        return 1;
    }
}
