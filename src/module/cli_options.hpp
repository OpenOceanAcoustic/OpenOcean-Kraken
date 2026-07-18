#ifndef OPENOCEANKRAKEN_CLI_OPTIONS_HPP
#define OPENOCEANKRAKEN_CLI_OPTIONS_HPP

#include <string>
#include <vector>

namespace OpenOceanKraken::cli
{
    struct Options
    {
        std::string input_path;
        std::string output_root;
        int threads = 0;
        bool velocity = false;
        bool export_mod = false;
        bool mod_only = false;
        bool export_json = false;
        bool help = false;
    };

    Options parse_options(const std::vector<std::string> &args);
    std::string usage();
}

#endif
