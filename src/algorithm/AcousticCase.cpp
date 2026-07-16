#include "algorithm/AcousticCase.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace OpenOceanKrakenc
{
namespace
{
std::string trim(std::string value)
{
    const auto notSpace = [](unsigned char c) { return !std::isspace(c); };
    value.erase(value.begin(), std::find_if(value.begin(), value.end(), notSpace));
    value.erase(std::find_if(value.rbegin(), value.rend(), notSpace).base(), value.end());
    return value;
}

std::string stripComment(const std::string &line)
{
    bool quoted = false;
    for (std::size_t i = 0; i < line.size(); ++i)
    {
        if (line[i] == '\'')
        {
            quoted = !quoted;
        }
        else if (line[i] == '!' && !quoted)
        {
            return line.substr(0, i);
        }
    }
    return line;
}

std::string quotedValue(const std::string &line, const char *field)
{
    const std::size_t first = line.find('\'');
    const std::size_t second = first == std::string::npos ? std::string::npos
                                                          : line.find('\'', first + 1);
    if (first == std::string::npos || second == std::string::npos)
    {
        throw std::runtime_error(std::string("missing quoted ") + field);
    }
    return line.substr(first + 1, second - first - 1);
}

std::vector<double> numbers(std::string line, const char *field)
{
    std::replace(line.begin(), line.end(), '/', ' ');
    std::replace(line.begin(), line.end(), ',', ' ');
    std::istringstream stream(line);
    std::vector<double> result;
    double value = 0.0;
    while (stream >> value)
    {
        result.push_back(value);
    }
    if (result.empty())
    {
        throw std::runtime_error(std::string("missing numeric ENV field: ") + field);
    }
    return result;
}

AcousticBoundaryType boundaryType(char value, const char *field)
{
    switch (value)
    {
    case 'V':
        return AcousticBoundaryType::Vacuum;
    case 'R':
        return AcousticBoundaryType::Rigid;
    case 'A':
        return AcousticBoundaryType::HalfSpace;
    case 'F':
        return AcousticBoundaryType::ReflectionCoefficient;
    case 'P':
        return AcousticBoundaryType::InternalReflection;
    default:
        throw std::runtime_error(std::string("unsupported phase-2 ") + field +
                                 " boundary: " + value);
    }
}

void validate(const AcousticCase &result)
{
    if (!(result.frequency > 0.0))
    {
        throw std::runtime_error("ENV frequency must be positive");
    }
    if (result.layers.empty())
    {
        throw std::runtime_error("ENV must contain an acoustic layer");
    }
    if (!(result.cLow > 0.0) || !(result.cHigh > result.cLow))
    {
        throw std::runtime_error("ENV phase-speed interval must satisfy 0 < cLow < cHigh");
    }
    bool hasAcousticLayer = false;
    for (const AcousticLayer &layer : result.layers)
    {
        if (layer.baseMesh < 1 || !(layer.bottomDepth > layer.topDepth) ||
            layer.samples.size() < 2)
        {
            throw std::runtime_error("invalid acoustic layer mesh or depth range");
        }
        if (std::abs(layer.samples.front().depth - layer.topDepth) > 1.0e-8 ||
            std::abs(layer.samples.back().depth - layer.bottomDepth) > 1.0e-8)
        {
            throw std::runtime_error("SSP samples must span the complete acoustic layer");
        }
        double previousDepth = layer.topDepth - 1.0;
        for (const AcousticSample &sample : layer.samples)
        {
            if (!(sample.depth > previousDepth) || !(sample.cp > 0.0) ||
                !(sample.rho > 0.0) || sample.alphaP < 0.0)
            {
                throw std::runtime_error("invalid or non-monotonic acoustic SSP sample");
            }
            if (sample.cs < 0.0 || sample.alphaS < 0.0)
            {
                throw std::runtime_error("invalid elastic SSP sample");
            }
            previousDepth = sample.depth;
        }
        hasAcousticLayer = hasAcousticLayer || layer.samples.front().cs == 0.0;
    }
    if (!hasAcousticLayer)
    {
        throw std::runtime_error("ENV must contain at least one acoustic layer");
    }
}

void readReflectionTable(const std::filesystem::path &path,
                         const char *extension,
                         AcousticBoundary &boundary)
{
    std::filesystem::path tablePath = path;
    tablePath.replace_extension(extension);
    std::ifstream table(tablePath);
    int count = 0;
    if (!(table >> count) || count < 2)
    {
        throw std::runtime_error("invalid or missing reflection file: " + tablePath.string());
    }
    boundary.reflectionSamples.clear();
    for (int index = 0; index < count; ++index)
    {
        ReflectionSample sample;
        double phaseDegrees = 0.0;
        if (!(table >> sample.angleDegrees >> sample.magnitude >> phaseDegrees))
        {
            throw std::runtime_error("truncated reflection file: " + tablePath.string());
        }
        sample.phaseRadians = phaseDegrees * std::acos(-1.0) / 180.0;
        boundary.reflectionSamples.push_back(sample);
    }
}
}

std::vector<AcousticCase> readAcousticEnvironments(
    const std::filesystem::path &path)
{
    std::ifstream input(path);
    if (!input)
    {
        throw std::runtime_error("unable to open ENV file: " + path.string());
    }

    std::vector<std::string> lines;
    std::string line;
    while (std::getline(input, line))
    {
        line = trim(stripComment(line));
        if (!line.empty())
        {
            lines.push_back(line);
        }
    }
    std::size_t cursor = 0;
    std::vector<AcousticCase> results;
    const auto take = [&](const char *field) -> const std::string & {
        if (cursor >= lines.size())
        {
            throw std::runtime_error(std::string("unexpected end of ENV while reading ") + field);
        }
        return lines[cursor++];
    };

    while (cursor < lines.size())
    {
    AcousticCase result;
    const auto readHalfSpace = [&](const char *field) {
        const std::vector<double> values = numbers(take(field), field);
        if (values.size() < 5)
        {
            throw std::runtime_error(std::string(field) +
                                     " requires depth, cp, cs, rho and alphaP");
        }
        AcousticBoundary boundary;
        boundary.type = AcousticBoundaryType::HalfSpace;
        boundary.depth = values[0];
        boundary.cp = values[1];
        boundary.cs = values[2];
        boundary.rho = values[3];
        boundary.alphaP = values[4];
        if (values.size() > 5)
        {
            boundary.alphaS = values[5];
        }
        return boundary;
    };
    result.title = quotedValue(take("title"), "title");
    result.frequency = numbers(take("frequency"), "frequency").front();
    result.referenceFrequency = result.frequency;

    const int mediumCount = static_cast<int>(numbers(take("medium count"), "medium count").front());
    if (mediumCount < 1)
    {
        throw std::runtime_error("ENV medium count must be positive");
    }

    const std::string options = quotedValue(take("top options"), "top options");
    if (options.size() < 3)
    {
        throw std::runtime_error("ENV top options must contain interpolation, boundary and attenuation");
    }
    switch (options[0])
    {
    case 'N':
        result.interpolation = AcousticInterpolation::N2Linear;
        break;
    case 'C':
        result.interpolation = AcousticInterpolation::CLinear;
        break;
    default:
        throw std::runtime_error("phase 2 supports only N2-linear and C-linear SSP interpolation");
    }
    result.top.type = boundaryType(options[1], "top");
    result.attenuationUnit = options[2];
    result.enableRootRestarts = options.size() >= 5 && options[4] == '.';
    if (result.top.type == AcousticBoundaryType::HalfSpace)
    {
        result.top = readHalfSpace("top half-space");
    }
    else if (result.top.type == AcousticBoundaryType::ReflectionCoefficient)
    {
        readReflectionTable(path, ".trc", result.top);
    }

    double topDepth = 0.0;
    for (int medium = 0; medium < mediumCount; ++medium)
    {
        const std::vector<double> header = numbers(take("medium header"), "medium header");
        if (header.size() < 3)
        {
            throw std::runtime_error("medium header requires mesh, roughness and bottom depth");
        }
        AcousticLayer layer;
        layer.baseMesh = static_cast<int>(header[0]);
        layer.topDepth = topDepth;
        layer.bottomDepth = header[2];

        AcousticSample previous;
        bool firstSample = true;
        while (layer.samples.empty() || layer.samples.back().depth < layer.bottomDepth)
        {
            const std::vector<double> values = numbers(take("SSP sample"), "SSP sample");
            if (values.size() < 2 || (firstSample && values.size() < 6))
            {
                throw std::runtime_error("first SSP sample requires depth, cp, cs, rho, alphaP and alphaS");
            }
            AcousticSample sample = previous;
            sample.depth = values[0];
            sample.cp = values[1];
            if (values.size() > 2)
                sample.cs = values[2];
            if (values.size() > 3)
                sample.rho = values[3];
            if (values.size() > 4)
                sample.alphaP = values[4];
            if (values.size() > 5)
                sample.alphaS = values[5];
            layer.samples.push_back(sample);
            previous = sample;
            firstSample = false;
        }
        result.layers.push_back(std::move(layer));
        topDepth = result.layers.back().bottomDepth;
    }

    const std::string bottomOptions = quotedValue(take("bottom options"), "bottom options");
    if (bottomOptions.empty())
    {
        throw std::runtime_error("missing bottom boundary option");
    }
    result.bottom.type = boundaryType(bottomOptions.front(), "bottom");
    result.bottom.depth = topDepth;
    if (result.bottom.type == AcousticBoundaryType::HalfSpace)
    {
        result.bottom = readHalfSpace("bottom half-space");
    }
    else if (result.bottom.type == AcousticBoundaryType::ReflectionCoefficient)
    {
        readReflectionTable(path, ".brc", result.bottom);
    }
    else if (result.bottom.type == AcousticBoundaryType::InternalReflection)
    {
        std::filesystem::path tablePath = path;
        tablePath.replace_extension(".irc");
        std::ifstream table(tablePath);
        std::string titleLine;
        int count = 0;
        if (!std::getline(table, titleLine) || !(table >> count) || count < 2)
        {
            throw std::runtime_error("invalid or missing IRC file: " + tablePath.string());
        }
        for (int index = 0; index < count; ++index)
        {
            InternalReflectionSample sample;
            double fReal = 0.0;
            double fImaginary = 0.0;
            double gReal = 0.0;
            double gImaginary = 0.0;
            if (!(table >> sample.eigenvalue >> fReal >> fImaginary >>
                  gReal >> gImaginary >> sample.power10))
            {
                throw std::runtime_error("truncated IRC file: " + tablePath.string());
            }
            sample.f = {fReal, fImaginary};
            sample.g = {gReal, gImaginary};
            result.bottom.internalSamples.push_back(sample);
        }
    }

    const std::vector<double> phaseSpeeds = numbers(take("phase-speed interval"), "phase-speed interval");
    if (phaseSpeeds.size() < 2)
    {
        throw std::runtime_error("phase-speed interval requires cLow and cHigh");
    }
    result.cLow = phaseSpeeds[0];
    result.cHigh = phaseSpeeds[1];
    result.rMaxKm = numbers(take("RMax"), "RMax").front();

    const auto readDepthVector = [&](const char *countField,
                                     const char *valuesField) {
        const int count = static_cast<int>(
            numbers(take(countField), countField).front());
        if (count < 1)
        {
            throw std::runtime_error(std::string(countField) + " must be positive");
        }
        const std::vector<double> values = numbers(take(valuesField), valuesField);
        if (static_cast<int>(values.size()) == count)
        {
            return values;
        }
        if (values.size() == 1)
        {
            return std::vector<double>(static_cast<std::size_t>(count), values.front());
        }
        if (values.size() == 2 && count > 1)
        {
            std::vector<double> expanded(static_cast<std::size_t>(count));
            for (int index = 0; index < count; ++index)
            {
                expanded[static_cast<std::size_t>(index)] =
                    values.front() + (values.back() - values.front()) *
                                         index / static_cast<double>(count - 1);
            }
            return expanded;
        }
        throw std::runtime_error(std::string(valuesField) + " count mismatch");
    };
    if (cursor < lines.size())
    {
        result.sourceDepths = readDepthVector("source depth count", "source depths");
    }
    if (cursor < lines.size())
    {
        result.receiverDepths = readDepthVector("receiver depth count", "receiver depths");
    }

    validate(result);
    results.push_back(std::move(result));
    }
    return results;
}

AcousticCase readAcousticEnv(const std::filesystem::path &path)
{
    std::vector<AcousticCase> profiles = readAcousticEnvironments(path);
    if (profiles.empty())
    {
        throw std::runtime_error("ENV file contains no profiles: " + path.string());
    }
    return std::move(profiles.front());
}

double soundSpeedAt(const AcousticCase &input, double depth)
{
    if (input.layers.empty())
    {
        throw std::invalid_argument("cannot interpolate sound speed in an empty profile");
    }
    if (depth <= input.layers.front().topDepth)
    {
        return input.top.cp > 0.0 ? input.top.cp
                                  : input.layers.front().samples.front().cp;
    }
    if (depth >= input.layers.back().bottomDepth)
    {
        return input.bottom.cp > 0.0 ? input.bottom.cp
                                     : input.layers.back().samples.back().cp;
    }
    for (const AcousticLayer &layer : input.layers)
    {
        if (depth > layer.bottomDepth)
        {
            continue;
        }
        const auto upper = std::upper_bound(
            layer.samples.begin(), layer.samples.end(), depth,
            [](double target, const AcousticSample &sample) {
                return target < sample.depth;
            });
        if (upper == layer.samples.begin())
        {
            return upper->cp;
        }
        if (upper == layer.samples.end())
        {
            return layer.samples.back().cp;
        }
        const AcousticSample &right = *upper;
        const AcousticSample &left = *(upper - 1);
        const double weight = (depth - left.depth) / (right.depth - left.depth);
        return left.cp + weight * (right.cp - left.cp);
    }
    throw std::runtime_error("sound-speed interpolation did not locate a layer");
}
}
