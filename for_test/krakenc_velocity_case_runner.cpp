#include "OpenOceanKrakencInterface.h"

#include <nlohmann/json.hpp>

#include <algorithm>
#include <cctype>
#include <exception>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>

namespace
{
int parseThreadCount(const char *value)
{
    const std::string text(value);
    if (text.empty() || !std::all_of(text.begin(), text.end(), [](unsigned char ch) {
            return std::isdigit(ch) != 0;
        }))
    {
        throw std::invalid_argument("thread count must be a positive integer");
    }
    const int threads = std::stoi(text);
    if (threads < 1 || threads > 1024)
    {
        throw std::invalid_argument("thread count must be between 1 and 1024");
    }
    return threads;
}
}

int main(int argc, char **argv)
{
    if (argc != 3 && argc != 4)
    {
        std::cerr << "Usage: OpenOceanKrakenc_velocity_case_runner "
                     "<case.env> <output-root> [threads]\n";
        return 2;
    }
    try
    {
        const std::filesystem::path envPath = argv[1];
        const std::filesystem::path outputRoot = argv[2];
        const int threadCount = argc == 4 ? parseThreadCount(argv[3]) : 1;
        if (envPath.extension() != ".env")
        {
            throw std::invalid_argument("velocity runner input must be an ENV file");
        }
        if (!outputRoot.parent_path().empty())
        {
            std::filesystem::create_directories(outputRoot.parent_path());
        }

        OpenOceanKrakenc::Interface api;
        const OpenOceanKrakenc::LoadResult loaded = api.loadEnv(envPath.string());
        if (!loaded.ok)
        {
            throw std::runtime_error(loaded.error.message);
        }
        api.setFieldThreads(threadCount);
        api.set_Velocity_enable(true);
        api.runEigen();
        api.runField();
        api.export_result(outputRoot.string());

        const OpenOceanKrakenc::OOKC_output &result = api.getOutput_const();
        nlohmann::ordered_json output;
        output["passes"] = true;
        output["source_count"] = result.sourceCount;
        output["receiver_depth_count"] = result.receiverDepthCount;
        output["range_count"] = result.rangeCount;
        output["threads"] = threadCount;
        output["pressure"] =
            std::filesystem::absolute(outputRoot.string() + "_P.shd").string();
        output["horizontal_velocity"] =
            std::filesystem::absolute(outputRoot.string() + "_H.shd").string();
        output["vertical_velocity"] =
            std::filesystem::absolute(outputRoot.string() + "_V.shd").string();
        output["mode_file"] =
            std::filesystem::absolute(outputRoot.string() + "_P.mod").string();
        std::cout << output.dump() << '\n';
        return 0;
    }
    catch (const std::exception &error)
    {
        std::cerr << "OpenOcean-Krakenc velocity export failed: "
                  << error.what() << '\n';
        return 1;
    }
}
