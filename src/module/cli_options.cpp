#include "cli_options.hpp"

#include <stdexcept>

namespace OpenOceanKraken::cli
{
    namespace
    {
        int parse_positive_integer(const std::string &value, const std::string &option)
        {
            std::size_t consumed = 0;
            int parsed = 0;
            try
            {
                parsed = std::stoi(value, &consumed);
            }
            catch (const std::exception &)
            {
                throw std::invalid_argument(option + " requires a positive integer");
            }
            if (consumed != value.size() || parsed <= 0)
            {
                throw std::invalid_argument(option + " requires a positive integer");
            }
            return parsed;
        }

        const std::string &require_value(
            const std::vector<std::string> &args,
            std::size_t &index,
            const std::string &option)
        {
            if (index + 1 >= args.size() || args[index + 1].rfind("--", 0) == 0)
            {
                throw std::invalid_argument(option + " requires a value");
            }
            return args[++index];
        }
    }

    Options parse_options(const std::vector<std::string> &args)
    {
        Options options;
        for (std::size_t index = 0; index < args.size(); ++index)
        {
            const std::string &argument = args[index];
            if (argument == "--help" || argument == "-h")
            {
                options.help = true;
            }
            else if (argument == "--output")
            {
                options.output_root = require_value(args, index, argument);
            }
            else if (argument == "--threads")
            {
                options.threads = parse_positive_integer(
                    require_value(args, index, argument), argument);
            }
            else if (argument == "--velocity")
            {
                options.velocity = true;
            }
            else if (argument == "--mod")
            {
                options.export_mod = true;
            }
            else if (argument == "--mod-only")
            {
                options.mod_only = true;
                options.export_mod = true;
            }
            else if (argument == "--json")
            {
                options.export_json = true;
            }
            else if (!argument.empty() && argument.front() == '-')
            {
                throw std::invalid_argument("Unknown option: " + argument);
            }
            else if (options.input_path.empty())
            {
                options.input_path = argument;
            }
            else
            {
                throw std::invalid_argument("Only one input file may be provided");
            }
        }
        return options;
    }

    std::string usage()
    {
        return
            "Usage: OpenOceanKraken <input.env|input.json> [options]\n"
            "Options:\n"
            "  --output <path>   Result root without an extension\n"
            "  --threads <N>     Positive worker-thread count\n"
            "  --velocity        Export pressure and vertical/horizontal velocity SHD\n"
            "  --mod             Export MOD in addition to the field result\n"
            "  --mod-only        Run eigenmode calculation and export only MOD\n"
            "  --json            Export the normalized input as JSON\n"
            "  --help, -h        Show this help text\n";
    }
}
