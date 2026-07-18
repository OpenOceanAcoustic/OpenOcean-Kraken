#include "OpenOceanKrakenInterface.h"
#include "ThreadPool.h"
#include "cli_options.hpp"

#include <algorithm>
#include <cctype>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

namespace
{
    std::string lower_extension(const std::filesystem::path &path)
    {
        std::string extension = path.extension().string();
        std::transform(
            extension.begin(),
            extension.end(),
            extension.begin(),
            [](unsigned char value) { return static_cast<char>(std::tolower(value)); });
        return extension;
    }

    std::string default_output_root(const std::filesystem::path &input)
    {
        std::filesystem::path root = input;
        root.replace_extension();
        return root.string();
    }

    int run_kraken(const OpenOceanKraken::cli::Options &options)
    {
        const std::filesystem::path input_path(options.input_path);
        if (!std::filesystem::is_regular_file(input_path))
        {
            throw std::runtime_error("Input file does not exist: " + input_path.string());
        }

        const std::string extension = lower_extension(input_path);
        if (extension != ".env" && extension != ".json")
        {
            throw std::invalid_argument("Input must have an .env or .json extension");
        }

        const std::string output_root = options.output_root.empty()
            ? default_output_root(input_path)
            : options.output_root;
        const std::filesystem::path output_path(output_root);
        if (!output_path.parent_path().empty())
        {
            std::filesystem::create_directories(output_path.parent_path());
        }

        const unsigned int hardware_threads = std::thread::hardware_concurrency();
        const int threads = options.threads > 0
            ? options.threads
            : static_cast<int>(std::max(1u, hardware_threads));

        ThreadPool thread_pool(static_cast<std::size_t>(threads));
        OpenOceanKraken::Interface interface(thread_pool);
        interface.setNumThreads(threads);

        const bool loaded = extension == ".json"
            ? interface.from_json(input_path.string())
            : interface.from_env(input_path.string());
        if (!loaded)
        {
            throw std::runtime_error("Failed to load input: " + input_path.string());
        }

        if (options.velocity)
        {
            interface.set_Velocity_enable(true);
        }
        if (options.export_json && !interface.to_json(output_root + ".json"))
        {
            throw std::runtime_error("Failed to export JSON: " + output_root + ".json");
        }

        if (options.mod_only)
        {
            interface.runEigen();
            interface.export_mod(output_root);
            std::cout << "OOK MOD written: " << output_root << ".mod\n";
            return 0;
        }

        interface.run();
        if (options.velocity)
        {
            interface.export_result(output_root);
            std::cout << "OOK pressure SHD written: " << output_root << "_P.shd\n";
            std::cout << "OOK vertical velocity SHD written: " << output_root << "_V.shd\n";
            std::cout << "OOK horizontal velocity SHD written: " << output_root << "_H.shd\n";
        }
        else
        {
            interface.export_shd(output_root, 1);
            std::cout << "OOK pressure SHD written: " << output_root << ".shd\n";
        }
        if (options.export_mod)
        {
            interface.export_mod(output_root);
            std::cout << "OOK MOD written: " << output_root << ".mod\n";
        }
        return 0;
    }
}

int main(int argc, char **argv)
{
    try
    {
        const std::vector<std::string> args(argv + 1, argv + argc);
        const auto options = OpenOceanKraken::cli::parse_options(args);
        if (options.help)
        {
            std::cout << OpenOceanKraken::cli::usage();
            return 0;
        }
        if (options.input_path.empty())
        {
            std::cerr << OpenOceanKraken::cli::usage();
            return 2;
        }
        return run_kraken(options);
    }
    catch (const std::invalid_argument &error)
    {
        std::cerr << "OpenOceanKraken: " << error.what() << '\n';
        return 2;
    }
    catch (const std::exception &error)
    {
        std::cerr << "OpenOceanKraken failed: " << error.what() << '\n';
        return 1;
    }
}
