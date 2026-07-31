#include "OpenOceanKrakenKernelInterface.h"
#include "ThreadPool.h"
#include "cli_options.hpp"

#include <algorithm>
#include <cctype>
#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

namespace
{
    struct MainOptions
    {
        OpenOceanKraken::cli::Options ook;
        std::string output_directory;
        bool print_timing = false;
    };

    const std::string &require_value(
        const std::vector<std::string> &args,
        std::size_t &index,
        const std::string &option)
    {
        if (index + 1 >= args.size() ||
            (!args[index + 1].empty() && args[index + 1].front() == '-'))
        {
            throw std::invalid_argument(option + " requires a value");
        }
        return args[++index];
    }

    MainOptions parse_main_options(const std::vector<std::string> &args)
    {
        MainOptions result;
        std::vector<std::string> normalized;
        normalized.reserve(args.size());

        for (std::size_t index = 0; index < args.size(); ++index)
        {
            const std::string &argument = args[index];
            if (argument == "-t")
            {
                normalized.emplace_back("--threads");
                normalized.push_back(require_value(args, index, argument));
            }
            else if (argument == "-v")
            {
                normalized.emplace_back("--velocity");
            }
            else if (argument == "-j")
            {
                normalized.emplace_back("--json");
            }
            else if (argument == "-out")
            {
                if (!result.output_directory.empty())
                {
                    throw std::invalid_argument("-out may only be specified once");
                }
                result.output_directory = require_value(args, index, argument);
            }
            else if (argument == "-time" || argument == "--time")
            {
                result.print_timing = true;
            }
            else if (argument == "-M")
            {
                throw std::invalid_argument(
                    "-M is not supported because OOK has no memory-limit interface");
            }
            else if (argument == "-m")
            {
                throw std::invalid_argument(
                    "-m is not supported because OOK has no memory-tracker interface");
            }
            else
            {
                normalized.push_back(argument);
            }
        }

        result.ook = OpenOceanKraken::cli::parse_options(normalized);
        if (!result.output_directory.empty() && !result.ook.output_root.empty())
        {
            throw std::invalid_argument("-out and --output cannot be used together");
        }
        return result;
    }

    std::string usage(const std::string &program_name)
    {
        return
            "Usage: " + program_name + " <input.env|input.json> [options]\n"
            "用法: " + program_name +
                " <文件名> [-t 线程数] [-v] [-j] [-out 输出目录] [-time]\n"
            "\n"
            "参数说明:\n"
            "  <文件名>            输入配置文件，支持 .env 或 .json；省略后缀时默认使用 .env\n"
            "  -t <线程数>         指定线程数量，等价于 --threads <N>\n"
            "  -v                  启用垂直、水平振速计算，等价于 --velocity\n"
            "  -j                  运行前导出 JSON 配置，等价于 --json\n"
            "  -out <输出目录>     将结果写入指定目录\n"
            "  -time, --time       打印计算耗时\n"
            "  --output <路径>     指定不含扩展名的结果根路径\n"
            "  --mod               在声场结果之外导出 MOD 文件\n"
            "  --mod-only          只计算本征模态并导出 MOD 文件\n"
            "  -h, --help          显示此帮助信息\n"
            "\n"
            "示例:\n"
            "  默认调用:           " + program_name + " MunkK.env\n"
            "  指定线程数:         " + program_name + " MunkK.env -t 8\n"
            "  启用振速计算:       " + program_name + " MunkK.env -v\n"
            "  导出 JSON:          " + program_name + " MunkK.env -j\n"
            "  指定输出目录:       " + program_name + " MunkK.env -out output\n"
            "  打印计算耗时:       " + program_name + " MunkK.env -time\n"
            "  组合调用:           " + program_name +
                " MunkK.env -t 8 -v -j -out output -time\n";
    }

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

    std::filesystem::path resolve_input_path(const std::string &input)
    {
        std::filesystem::path input_path(input);
        if (input_path.extension().empty())
        {
            input_path += ".env";
        }

        const std::string extension = lower_extension(input_path);
        if (extension != ".env" && extension != ".json")
        {
            throw std::invalid_argument("Input must have an .env or .json extension");
        }
        if (!std::filesystem::is_regular_file(input_path))
        {
            throw std::runtime_error("Input file does not exist: " + input_path.string());
        }
        return input_path;
    }

    std::string default_output_root(const std::filesystem::path &input)
    {
        std::filesystem::path root = input;
        root.replace_extension();
        return root.string();
    }

    std::string select_output_root(
        const MainOptions &options,
        const std::filesystem::path &input_path)
    {
        if (!options.ook.output_root.empty())
        {
            return options.ook.output_root;
        }
        if (!options.output_directory.empty())
        {
            return (
                std::filesystem::path(options.output_directory) /
                input_path.stem()).string();
        }
        return default_output_root(input_path);
    }

    template <typename Function>
    void run_with_optional_timing(bool print_timing, Function &&function)
    {
        if (!print_timing)
        {
            function();
            return;
        }

        const auto start = std::chrono::steady_clock::now();
        function();
        const auto finish = std::chrono::steady_clock::now();
        const auto milliseconds =
            std::chrono::duration_cast<std::chrono::milliseconds>(finish - start).count();
        std::cout << "OOK calculation time: " << milliseconds << " milliseconds\n";
    }

    int run_kraken(const MainOptions &options)
    {
        const std::filesystem::path input_path =
            resolve_input_path(options.ook.input_path);
        const std::string extension = lower_extension(input_path);
        const std::string output_root = select_output_root(options, input_path);
        const std::filesystem::path output_path(output_root);
        if (!output_path.parent_path().empty())
        {
            std::filesystem::create_directories(output_path.parent_path());
        }

        const unsigned int hardware_threads = std::thread::hardware_concurrency();
        const int threads = options.ook.threads > 0
            ? options.ook.threads
            : static_cast<int>(std::max(1u, hardware_threads));

        ThreadPool thread_pool(static_cast<std::size_t>(threads));
        OpenOceanKraken::KernelInterface interface(thread_pool);
        interface.setNumThreads(threads);

        const bool loaded = extension == ".json"
            ? interface.from_json(input_path.string())
            : interface.from_env(input_path.string());
        if (!loaded)
        {
            throw std::runtime_error("Failed to load input: " + input_path.string());
        }

        if (options.ook.velocity)
        {
            interface.set_Velocity_enable(true);
        }
        if (options.ook.export_json && !interface.to_json(output_root + ".json"))
        {
            throw std::runtime_error("Failed to export JSON: " + output_root + ".json");
        }

        std::cout << "numThreads: " << threads << '\n';

        if (options.ook.mod_only)
        {
            run_with_optional_timing(
                options.print_timing,
                [&interface]() { interface.runEigen(); });
            interface.export_mod(output_root);
            std::cout << "OOK MOD written: " << output_root << ".mod\n";
            std::cout << "OpenOcean-Kraken done.\n";
            return 0;
        }

        run_with_optional_timing(
            options.print_timing,
            [&interface]() { interface.run(); });
        if (options.ook.velocity)
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
        if (options.ook.export_mod)
        {
            interface.export_mod(output_root);
            std::cout << "OOK MOD written: " << output_root << ".mod\n";
        }
        std::cout << "OpenOcean-Kraken done.\n";
        return 0;
    }
}

int main(int argc, char **argv)
{
#ifdef _WIN32
    std::system("chcp 65001 >nul");
#endif

    try
    {
        const std::string program_name =
            argc > 0 && argv[0] != nullptr ? argv[0] : "OpenOceanKraken";
        if (argc == 1)
        {
            std::cerr << usage(program_name);
            return 2;
        }

        // 与 OOB 一致：第一个参数为 -h 或 --help 时直接显示完整帮助。
        const std::string first_argument = argv[1];
        if (first_argument == "-h" || first_argument == "--help")
        {
            std::cout << usage(program_name);
            return 0;
        }

        const std::vector<std::string> args(argv + 1, argv + argc);
        const MainOptions options = parse_main_options(args);
        if (options.ook.help)
        {
            std::cout << usage(program_name);
            return 0;
        }
        if (options.ook.input_path.empty())
        {
            std::cerr << usage(program_name);
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
