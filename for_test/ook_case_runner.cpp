#include "OpenOceanKrakenInterface.h"
#include "ThreadPool.h"

#include <algorithm>
#include <exception>
#include <iostream>
#include <string>
#include <thread>

int main(int argc, char **argv)
{
    if (argc < 3)
    {
        std::cerr << "Usage: ook_case_runner <env_path> <output_root_without_extension> [threads] [--mod]\n";
        return 2;
    }

    const std::string env_path = argv[1];
    const std::string output_root = argv[2];
    int threads = 1;
    bool export_mod = false;
    if (argc >= 4)
    {
        for (int i = 3; i < argc; ++i)
        {
            const std::string arg = argv[i];
            if (arg == "--mod")
            {
                export_mod = true;
            }
            else
            {
                threads = std::max(1, std::stoi(arg));
            }
        }
    }
    else
    {
        const unsigned int hw = std::thread::hardware_concurrency();
        threads = std::max(1u, hw == 0 ? 1u : hw);
    }

    try
    {
        ThreadPool thread_pool(static_cast<size_t>(threads));
        OpenOceanKraken::Interface kraken_interface(thread_pool);
        kraken_interface.setNumThreads(threads);

        if (!kraken_interface.from_env(env_path))
        {
            std::cerr << "Failed to load ENV file: " << env_path << "\n";
            return 1;
        }

        kraken_interface.run();
        if (export_mod)
        {
            kraken_interface.export_mod(output_root);
            std::cout << "OOK MOD written: " << output_root << ".mod\n";
        }
        kraken_interface.export_shd(output_root, 1);
        std::cout << "OOK SHD written: " << output_root << ".shd\n";
        return 0;
    }
    catch (const std::exception &e)
    {
        std::cerr << "OOK runner failed: " << e.what() << "\n";
        return 1;
    }
}
