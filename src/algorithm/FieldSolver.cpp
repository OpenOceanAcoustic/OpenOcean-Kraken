#include "algorithm/FieldSolver.h"

#include "algorithm/ComplexNumerics.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstring>
#include <fstream>
#include <future>
#include <iterator>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace OpenOceanKrakenc
{
namespace
{
template <typename Value>
Value getValue(const std::vector<char> &data, std::size_t offset)
{
    if (offset + sizeof(Value) > data.size())
    {
        throw std::runtime_error("binary acoustic file is truncated");
    }
    Value result;
    std::memcpy(&result, data.data() + offset, sizeof(Value));
    return result;
}

template <typename Value>
void putValue(std::vector<char> &record, std::size_t offset, const Value &value)
{
    if (offset + sizeof(Value) > record.size())
    {
        throw std::runtime_error("SHD record overflow");
    }
    std::memcpy(record.data() + offset, &value, sizeof(Value));
}

void writeRecord(std::ofstream &stream, std::size_t index,
                 const std::vector<char> &record)
{
    stream.seekp(static_cast<std::streamoff>(index * record.size()));
    stream.write(record.data(), static_cast<std::streamsize>(record.size()));
    if (!stream)
    {
        throw std::runtime_error("failed while writing SHD record");
    }
}

std::string trim(std::string value)
{
    const auto nonSpace = [](unsigned char c) { return !std::isspace(c); };
    value.erase(value.begin(), std::find_if(value.begin(), value.end(), nonSpace));
    value.erase(std::find_if(value.rbegin(), value.rend(), nonSpace).base(), value.end());
    return value;
}

std::string cleanLine(const std::string &line)
{
    bool quoted = false;
    for (std::size_t index = 0; index < line.size(); ++index)
    {
        if (line[index] == '\'')
        {
            quoted = !quoted;
        }
        else if (line[index] == '!' && !quoted)
        {
            return trim(line.substr(0, index));
        }
    }
    return trim(line);
}

std::string unquote(std::string value)
{
    value = trim(value);
    const std::size_t first = value.find('\'');
    const std::size_t last = value.rfind('\'');
    if (first != std::string::npos && last > first)
    {
        return value.substr(first + 1, last - first - 1);
    }
    if (!value.empty() && value.front() == '/')
    {
        return {};
    }
    return value;
}

std::vector<double> numbers(std::string line)
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
    return result;
}

std::vector<double> expandVector(int count, const std::vector<double> &values)
{
    if (count < 1 || values.empty())
    {
        throw std::runtime_error("FLP vector is empty");
    }
    if (static_cast<int>(values.size()) == count)
    {
        return values;
    }
    if (values.size() == 1)
    {
        return std::vector<double>(static_cast<std::size_t>(count), values.front());
    }
    if (values.size() == 2)
    {
        if (count == 1)
        {
            return {values.front()};
        }
        std::vector<double> result(static_cast<std::size_t>(count));
        for (int index = 0; index < count; ++index)
        {
            result[static_cast<std::size_t>(index)] =
                values.front() + (values.back() - values.front()) *
                                     index / static_cast<double>(count - 1);
        }
        return result;
    }
    throw std::runtime_error("FLP vector count mismatch");
}

std::complex<double> interpolate(const ModeProfileData &data,
                                 std::size_t mode,
                                 double target)
{
    if (target <= data.depths.front())
    {
        return data.modes[mode].front();
    }
    if (target >= data.depths.back())
    {
        return data.modes[mode].back();
    }
    const auto right = std::upper_bound(data.depths.begin(), data.depths.end(), target);
    const std::size_t rightIndex = static_cast<std::size_t>(right - data.depths.begin());
    const std::size_t leftIndex = rightIndex - 1;
    const double fraction = (target - data.depths[leftIndex]) /
                            (data.depths[rightIndex] - data.depths[leftIndex]);
    return (1.0 - fraction) * data.modes[mode][leftIndex] +
           fraction * data.modes[mode][rightIndex];
}

std::complex<double> roundedComplex(std::complex<double> value)
{
    const std::complex<float> rounded(static_cast<float>(value.real()),
                                      static_cast<float>(value.imag()));
    return {rounded.real(), rounded.imag()};
}

std::complex<double> interpolateModeAtDepth(
    const ModeProfileData &profile,
    std::size_t mode,
    double target,
    double frequency)
{
    const auto halfspaceValue = [&](const ModeBoundaryData &boundary,
                                    double distance,
                                    std::complex<double> boundaryMode) {
        if (boundary.type != 'A' || boundary.cp == std::complex<double>(0.0, 0.0))
        {
            return boundaryMode;
        }
        constexpr float piMode = 3.141592f;
        const std::complex<double> mediumK =
            static_cast<double>(2.0f * piMode * static_cast<float>(frequency)) /
            boundary.cp;
        const std::complex<double> gamma = roundedComplex(pekerisRoot(
            profile.wavenumbers[mode] * profile.wavenumbers[mode] -
            mediumK * mediumK));
        return roundedComplex(boundaryMode * std::exp(-gamma * distance));
    };
    if (target < profile.top.depth)
    {
        return halfspaceValue(profile.top, profile.top.depth - target,
                              profile.modes[mode].front());
    }
    if (target > profile.bottom.depth)
    {
        return halfspaceValue(profile.bottom, target - profile.bottom.depth,
                              profile.modes[mode].back());
    }
    return roundedComplex(interpolate(profile, mode, target));
}

std::complex<double> interpolateModeDerivativeAtDepth(
    const ModeProfileData &profile,
    std::size_t mode,
    double target,
    double frequency)
{
    if (target < profile.top.depth && profile.top.type == 'A')
    {
        const std::complex<double> value = interpolateModeAtDepth(
            profile, mode, target, frequency);
        const std::complex<double> mediumK =
            (2.0 * pi * frequency) / profile.top.cp;
        return pekerisRoot(profile.wavenumbers[mode] * profile.wavenumbers[mode] -
                           mediumK * mediumK) * value;
    }
    if (target > profile.bottom.depth && profile.bottom.type == 'A')
    {
        const std::complex<double> value = interpolateModeAtDepth(
            profile, mode, target, frequency);
        const std::complex<double> mediumK =
            (2.0 * pi * frequency) / profile.bottom.cp;
        return -pekerisRoot(profile.wavenumbers[mode] * profile.wavenumbers[mode] -
                            mediumK * mediumK) * value;
    }
    if (profile.depths.size() < 2)
    {
        return 0.0;
    }
    std::size_t right = 1;
    if (target >= profile.depths.back())
    {
        right = profile.depths.size() - 1;
    }
    else if (target > profile.depths.front())
    {
        right = static_cast<std::size_t>(std::upper_bound(
            profile.depths.begin(), profile.depths.end(), target) -
            profile.depths.begin());
    }
    const std::size_t left = right - 1;
    return (profile.modes[mode][right] - profile.modes[mode][left]) /
           (profile.depths[right] - profile.depths[left]);
}

double interpolateBeam(const FieldParameters &parameters, double angle)
{
    if (parameters.beamAnglesDegrees.empty())
    {
        return 1.0;
    }
    if (angle <= parameters.beamAnglesDegrees.front())
    {
        return parameters.beamAmplitudes.front();
    }
    if (angle >= parameters.beamAnglesDegrees.back())
    {
        return parameters.beamAmplitudes.back();
    }
    const auto right = std::upper_bound(
        parameters.beamAnglesDegrees.begin(), parameters.beamAnglesDegrees.end(), angle);
    const std::size_t rightIndex =
        static_cast<std::size_t>(right - parameters.beamAnglesDegrees.begin());
    const std::size_t leftIndex = rightIndex - 1;
    const double fraction = (angle - parameters.beamAnglesDegrees[leftIndex]) /
                            (parameters.beamAnglesDegrees[rightIndex] -
                             parameters.beamAnglesDegrees[leftIndex]);
    return (1.0 - fraction) * parameters.beamAmplitudes[leftIndex] +
           fraction * parameters.beamAmplitudes[rightIndex];
}

template <typename Function>
void parallelFor(std::size_t count, int requestedThreads, Function function)
{
    const std::size_t workers = std::min(
        count, static_cast<std::size_t>(std::max(1, requestedThreads)));
    if (workers <= 1)
    {
        for (std::size_t index = 0; index < count; ++index)
        {
            function(index);
        }
        return;
    }
    for (std::size_t batch = 0; batch < count; batch += workers)
    {
        std::vector<std::future<void>> futures;
        const std::size_t end = std::min(count, batch + workers);
        for (std::size_t index = batch; index < end; ++index)
        {
            futures.push_back(std::async(std::launch::async, function, index));
        }
        for (std::future<void> &future : futures)
        {
            future.get();
        }
    }
}
}

void validateFieldParameters(const FieldParameters &parameters)
{
    const std::string sourceTypes = "RXS";
    if (sourceTypes.find(parameters.sourceType) == std::string::npos)
    {
        throw std::invalid_argument("FLP source type must be R, X, or S");
    }
    if (parameters.propagationType != ' ' &&
        parameters.propagationType != 'A' &&
        parameters.propagationType != 'C')
    {
        throw std::invalid_argument("FLP propagation type must be A or C");
    }
    if (parameters.modeLimit < 1 || parameters.threadCount < 1)
    {
        throw std::invalid_argument("FLP mode limit and thread count must be positive");
    }
    if (parameters.profileRangesKm.empty() || parameters.rangesMetres.empty() ||
        parameters.sourceDepths.empty() || parameters.receiverDepths.empty() ||
        parameters.rangeOffsets.size() != parameters.receiverDepths.size())
    {
        throw std::invalid_argument("FLP range and depth vectors must be nonempty and consistent");
    }
    if (parameters.profileRangesKm.front() != 0.0)
    {
        throw std::invalid_argument("FLP first profile range must be zero");
    }
    for (std::size_t index = 0; index < parameters.profileRangesKm.size(); ++index)
    {
        if (!std::isfinite(parameters.profileRangesKm[index]) ||
            parameters.profileRangesKm[index] < 0.0 ||
            (index > 0 && parameters.profileRangesKm[index] <=
                              parameters.profileRangesKm[index - 1]))
        {
            throw std::invalid_argument("FLP profile ranges must be finite and strictly increasing");
        }
    }
    for (std::size_t index = 0; index < parameters.rangesMetres.size(); ++index)
    {
        if (!std::isfinite(parameters.rangesMetres[index]) ||
            parameters.rangesMetres[index] < 0.0 ||
            (index > 0 && parameters.rangesMetres[index] <
                              parameters.rangesMetres[index - 1]))
        {
            throw std::invalid_argument("FLP receiver ranges must be finite, nonnegative, and nondecreasing");
        }
    }
    const auto finiteValues = [](const std::vector<double> &values) {
        return std::all_of(values.begin(), values.end(), [](double value) {
            return std::isfinite(value);
        });
    };
    if (!finiteValues(parameters.sourceDepths) ||
        !finiteValues(parameters.receiverDepths) ||
        !finiteValues(parameters.rangeOffsets))
    {
        throw std::invalid_argument("FLP depths and receiver offsets must be finite");
    }
    if (!parameters.sourceSoundSpeeds.empty() &&
        (parameters.sourceSoundSpeeds.size() != parameters.sourceDepths.size() ||
         !std::all_of(parameters.sourceSoundSpeeds.begin(),
                      parameters.sourceSoundSpeeds.end(), [](double value) {
                          return std::isfinite(value) && value > 0.0;
                      })))
    {
        throw std::invalid_argument("source sound speeds must match source depths and be positive");
    }
    if (parameters.sourceType == 'R')
    {
        for (double offset : parameters.rangeOffsets)
        {
            if (parameters.rangesMetres.front() + offset < 0.0)
            {
                throw std::invalid_argument("cylindrical effective receiver ranges must be nonnegative");
            }
        }
    }
    if (parameters.beamPattern)
    {
        if (parameters.beamAnglesDegrees.size() < 2 ||
            parameters.beamAnglesDegrees.size() != parameters.beamAmplitudes.size() ||
            !finiteValues(parameters.beamAnglesDegrees) ||
            !finiteValues(parameters.beamAmplitudes))
        {
            throw std::invalid_argument("SBP angle and amplitude vectors are inconsistent");
        }
        for (std::size_t index = 1; index < parameters.beamAnglesDegrees.size(); ++index)
        {
            if (parameters.beamAnglesDegrees[index] <=
                parameters.beamAnglesDegrees[index - 1])
            {
                throw std::invalid_argument("SBP angles must be strictly increasing");
            }
        }
    }
}

FieldParameters readFieldParameters(const std::filesystem::path &path)
{
    std::ifstream stream(path);
    if (!stream)
    {
        throw std::runtime_error("unable to open FLP file: " + path.string());
    }
    std::vector<std::string> lines;
    std::string line;
    while (std::getline(stream, line))
    {
        line = cleanLine(line);
        if (!line.empty())
        {
            lines.push_back(line);
        }
    }
    std::size_t cursor = 0;
    const auto take = [&]() -> const std::string & {
        if (cursor >= lines.size())
        {
            throw std::runtime_error("unexpected end of FLP file");
        }
        return lines[cursor++];
    };
    const auto scalarInt = [&]() {
        const std::vector<double> value = numbers(take());
        if (value.empty())
        {
            throw std::runtime_error("missing FLP integer");
        }
        return static_cast<int>(value.front());
    };
    const auto vector = [&]() {
        const int count = scalarInt();
        return expandVector(count, numbers(take()));
    };

    FieldParameters result;
    result.title = unquote(take());
    const std::string options = unquote(take());
    if (options.empty())
    {
        throw std::runtime_error("missing FLP options");
    }
    result.sourceType = options[0];
    result.propagationType = options.size() >= 2 ? options[1] : ' ';
    if (options.size() >= 3 && options[2] != '*' && options[2] != 'O' &&
        options[2] != ' ')
    {
        throw std::runtime_error("FLP beam option must be '*', O, or blank");
    }
    if (options.size() >= 4 && options[3] != 'I' && options[3] != 'C' &&
        options[3] != ' ')
    {
        throw std::runtime_error("FLP coherence option must be I, C, or blank");
    }
    result.beamPattern = options.size() >= 3 && options[2] == '*';
    result.coherent = !(options.size() >= 4 && options[3] == 'I');
    result.modeLimit = scalarInt();
    result.profileRangesKm = vector();
    result.rangesMetres = vector();
    for (double &range : result.rangesMetres)
    {
        range *= 1000.0;
    }
    result.sourceDepths = vector();
    result.receiverDepths = vector();
    result.rangeOffsets = vector();
    if (result.rangeOffsets.size() != result.receiverDepths.size())
    {
        throw std::runtime_error("FLP receiver offsets must match receiver depths");
    }
    if (result.beamPattern)
    {
        std::filesystem::path beamPath = path;
        beamPath.replace_extension(".sbp");
        std::ifstream beam(beamPath);
        int count = 0;
        if (!(beam >> count) || count < 2)
        {
            throw std::runtime_error("invalid or missing SBP file: " + beamPath.string());
        }
        for (int index = 0; index < count; ++index)
        {
            double angle = 0.0;
            double decibels = 0.0;
            if (!(beam >> angle >> decibels))
            {
                throw std::runtime_error("truncated SBP file");
            }
            result.beamAnglesDegrees.push_back(angle);
            result.beamAmplitudes.push_back(std::pow(10.0, decibels / 20.0));
        }
    }
    validateFieldParameters(result);
    return result;
}

ModeFileData readModeFile(const std::filesystem::path &path)
{
    std::ifstream stream(path, std::ios::binary);
    if (!stream)
    {
        throw std::runtime_error("unable to open MOD file: " + path.string());
    }
    const std::vector<char> data(
        (std::istreambuf_iterator<char>(stream)), std::istreambuf_iterator<char>());
    const int recordWords = getValue<int>(data, 0);
    if (recordWords < 25)
    {
        throw std::runtime_error("invalid MOD record length");
    }
    const std::size_t recordBytes = static_cast<std::size_t>(4 * recordWords);
    ModeFileData result;
    std::size_t profileRecord = 0;
    bool firstProfile = true;
    while (profileRecord * recordBytes < data.size())
    {
        const std::size_t header = profileRecord * recordBytes;
        if (header + 100 > data.size())
        {
            throw std::runtime_error("MOD file is truncated in a profile header");
        }
        if (getValue<int>(data, header) != recordWords)
        {
            throw std::runtime_error("MOD profile record length mismatch");
        }
        const int frequencyCount = getValue<int>(data, header + 84);
        const int mediumCount = getValue<int>(data, header + 88);
        const int depthCount = getValue<int>(data, header + 92);
        const int materialCount = getValue<int>(data, header + 96);
        if (frequencyCount != 1 || mediumCount < 1 || depthCount < 1 ||
            materialCount != depthCount)
        {
            throw std::runtime_error(
                "MOD reader requires one-frequency acoustic mode profiles");
        }

        ModeProfileData profile;
        profile.meshCounts.resize(static_cast<std::size_t>(mediumCount));
        profile.materials.resize(static_cast<std::size_t>(mediumCount));
        profile.mediumDepths.resize(static_cast<std::size_t>(mediumCount));
        profile.mediumDensities.resize(static_cast<std::size_t>(mediumCount));
        for (int medium = 0; medium < mediumCount; ++medium)
        {
            const std::size_t materialOffset =
                (profileRecord + 1) * recordBytes + static_cast<std::size_t>(12 * medium);
            profile.meshCounts[static_cast<std::size_t>(medium)] =
                getValue<int>(data, materialOffset);
            profile.materials[static_cast<std::size_t>(medium)] = trim(
                std::string(data.data() + materialOffset + 4,
                            data.data() + materialOffset + 12));
            const std::size_t mediumOffset =
                (profileRecord + 2) * recordBytes + static_cast<std::size_t>(8 * medium);
            profile.mediumDepths[static_cast<std::size_t>(medium)] =
                getValue<float>(data, mediumOffset);
            profile.mediumDensities[static_cast<std::size_t>(medium)] =
                getValue<float>(data, mediumOffset + 4);
        }
        const double frequency = getValue<double>(data, (profileRecord + 3) * recordBytes);
        profile.depths.resize(static_cast<std::size_t>(depthCount));
        for (int index = 0; index < depthCount; ++index)
        {
            profile.depths[static_cast<std::size_t>(index)] = getValue<float>(
                data, (profileRecord + 4) * recordBytes + static_cast<std::size_t>(4 * index));
        }
        const int modeCount = getValue<int>(data, (profileRecord + 5) * recordBytes);
        if (modeCount <= 0)
        {
            throw std::runtime_error("MOD profile must contain at least one mode");
        }
        const auto readBoundary = [&](std::size_t &cursor) {
            ModeBoundaryData boundary;
            boundary.type = getValue<char>(data, cursor++);
            boundary.cp = {getValue<float>(data, cursor),
                           getValue<float>(data, cursor + 4)};
            cursor += 8;
            boundary.cs = {getValue<float>(data, cursor),
                           getValue<float>(data, cursor + 4)};
            cursor += 8;
            boundary.rho = getValue<float>(data, cursor);
            cursor += 4;
            boundary.depth = getValue<float>(data, cursor);
            cursor += 4;
            return boundary;
        };
        std::size_t boundaryCursor = (profileRecord + 6) * recordBytes;
        profile.top = readBoundary(boundaryCursor);
        profile.bottom = readBoundary(boundaryCursor);

        profile.modes.resize(static_cast<std::size_t>(modeCount));
        for (int mode = 0; mode < modeCount; ++mode)
        {
            profile.modes[static_cast<std::size_t>(mode)].resize(
                static_cast<std::size_t>(depthCount));
            const std::size_t start =
                (profileRecord + 7 + static_cast<std::size_t>(mode)) * recordBytes;
            for (int depth = 0; depth < depthCount; ++depth)
            {
                profile.modes[static_cast<std::size_t>(mode)]
                             [static_cast<std::size_t>(depth)] = {
                    getValue<float>(data, start + static_cast<std::size_t>(8 * depth)),
                    getValue<float>(data, start + static_cast<std::size_t>(8 * depth + 4))};
            }
        }
        profile.wavenumbers.resize(static_cast<std::size_t>(modeCount));
        const int complexPerRecord = recordWords / 2;
        const std::size_t waveStart = profileRecord + 7 + static_cast<std::size_t>(modeCount);
        for (int mode = 0; mode < modeCount; ++mode)
        {
            const std::size_t waveRecord =
                waveStart + static_cast<std::size_t>(mode / complexPerRecord);
            const std::size_t waveOffset = waveRecord * recordBytes +
                static_cast<std::size_t>(8 * (mode % complexPerRecord));
            profile.wavenumbers[static_cast<std::size_t>(mode)] = {
                getValue<float>(data, waveOffset), getValue<float>(data, waveOffset + 4)};
        }
        const std::size_t waveRecords = modeCount == 0 ? 0 :
            static_cast<std::size_t>((modeCount + complexPerRecord - 1) / complexPerRecord);
        profileRecord += 7 + static_cast<std::size_t>(modeCount) + waveRecords;

        if (firstProfile)
        {
            result.title.assign(data.data() + header + 4, data.data() + header + 84);
            result.title = trim(result.title);
            result.frequency = frequency;
            static_cast<ModeProfileData &>(result) = std::move(profile);
            firstProfile = false;
        }
        else
        {
            if (std::abs(frequency - result.frequency) > 1.0e-9)
            {
                throw std::runtime_error("multi-profile MOD frequencies do not match");
            }
            result.additionalProfiles.push_back(std::move(profile));
        }
    }
    if (firstProfile)
    {
        throw std::runtime_error("MOD file contains no profiles");
    }
    return result;
}

PressureField evaluateRangeIndependentField(
    const ModeFileData &modes,
    const FieldParameters &parameters)
{
    validateFieldParameters(parameters);
    if (modes.wavenumbers.size() != modes.modes.size() || modes.depths.empty() ||
        modes.modes.empty() || parameters.profileRangesKm.size() != 1)
    {
        throw std::invalid_argument("range-independent field inputs are inconsistent");
    }
    const std::size_t modeCount = std::min(
        modes.modes.size(), static_cast<std::size_t>(std::max(0, parameters.modeLimit)));
    PressureField result;
    result.title = parameters.title.empty() ? modes.title : parameters.title;
    result.frequency = modes.frequency;
    result.sourceDepths = parameters.sourceDepths;
    result.receiverDepths = parameters.receiverDepths;
    result.rangesMetres = parameters.rangesMetres;
    result.values.resize(result.sourceDepths.size() * result.receiverDepths.size() *
                         result.rangesMetres.size());
    result.horizontalValues.resize(result.values.size());
    result.verticalValues.resize(result.values.size());
    result.workerThreadsUsed = static_cast<int>(std::min(
        result.sourceDepths.size(),
        static_cast<std::size_t>(std::max(1, parameters.threadCount))));
    const std::complex<double> imaginary(0.0, 1.0);
    const std::complex<double> factor = imaginary * std::sqrt(2.0 * pi) *
                                        std::exp(imaginary * pi / 4.0);
    const double omega = 2.0 * pi * modes.frequency;
    std::vector<std::complex<double>> exponents(modeCount);
    std::vector<std::complex<double>> hank(
        result.rangesMetres.size() * modeCount);
    for (std::size_t mode = 0; mode < modeCount; ++mode)
    {
        exponents[mode] = -imaginary * modes.wavenumbers[mode];
        if (!parameters.coherent)
        {
            exponents[mode] = exponents[mode].real();
        }
        for (std::size_t range = 0; range < result.rangesMetres.size(); ++range)
        {
            hank[range * modeCount + mode] =
                std::exp(exponents[mode] * result.rangesMetres[range]);
        }
    }

    parallelFor(result.sourceDepths.size(), parameters.threadCount,
                [&](std::size_t source) {
        std::vector<std::complex<double>> coefficients(modeCount);
        std::vector<std::complex<double>> horizontalCoefficients(modeCount);
        std::vector<std::complex<double>> verticalCoefficients(modeCount);
        const double sourceSoundSpeed =
            source < parameters.sourceSoundSpeeds.size()
                ? parameters.sourceSoundSpeeds[source] : 1500.0;
        for (std::size_t mode = 0; mode < modeCount; ++mode)
        {
            std::complex<double> sourceMode = interpolateModeAtDepth(
                modes, mode, result.sourceDepths[source], modes.frequency);
            if (parameters.beamPattern && source == 0)
            {
                const double kz2 = std::max(
                    0.0, (omega * omega / (1500.0 * 1500.0) -
                          modes.wavenumbers[mode] * modes.wavenumbers[mode])
                             .real());
                const double angle = std::atan(
                    std::sqrt(kz2) / modes.wavenumbers[mode].real()) * 180.0 / pi;
                sourceMode *= interpolateBeam(parameters, angle);
            }
            if (parameters.sourceType == 'X')
            {
                coefficients[mode] = factor * sourceMode / modes.wavenumbers[mode];
                horizontalCoefficients[mode] = factor * sourceMode *
                    modes.wavenumbers[mode] * sourceSoundSpeed / omega;
                verticalCoefficients[mode] = factor * sourceMode * sourceSoundSpeed /
                    (modes.wavenumbers[mode] * omega * imaginary);
            }
            else
            {
                const std::complex<double> root = std::sqrt(modes.wavenumbers[mode]);
                coefficients[mode] = factor * sourceMode / root;
                horizontalCoefficients[mode] = factor * sourceMode * root *
                    sourceSoundSpeed / omega;
                verticalCoefficients[mode] = factor * sourceMode * sourceSoundSpeed /
                    (root * omega * imaginary);
            }
        }

        for (std::size_t receiver = 0; receiver < result.receiverDepths.size(); ++receiver)
        {
            std::vector<std::complex<double>> receiverCoefficients(modeCount);
            std::vector<std::complex<double>> receiverHorizontal(modeCount);
            std::vector<std::complex<double>> receiverVertical(modeCount);
            for (std::size_t mode = 0; mode < modeCount; ++mode)
            {
                const std::complex<double> offset =
                    std::exp(exponents[mode] * parameters.rangeOffsets[receiver]);
                const std::complex<double> receiverMode =
                    interpolateModeAtDepth(modes, mode,
                                           result.receiverDepths[receiver],
                                           modes.frequency);
                receiverCoefficients[mode] =
                    coefficients[mode] * receiverMode * offset;
                receiverHorizontal[mode] =
                    horizontalCoefficients[mode] * receiverMode * offset;
                receiverVertical[mode] = verticalCoefficients[mode] *
                    interpolateModeDerivativeAtDepth(
                        modes, mode, result.receiverDepths[receiver],
                        modes.frequency) *
                    offset;
            }
            for (std::size_t range = 0; range < result.rangesMetres.size(); ++range)
            {
                std::complex<double> pressure = 0.0;
                std::complex<double> horizontal = 0.0;
                std::complex<double> vertical = 0.0;
                for (std::size_t mode = 0; mode < modeCount; ++mode)
                {
                    const std::complex<double> propagation =
                        hank[range * modeCount + mode];
                    const std::complex<double> term = receiverCoefficients[mode] * propagation;
                    const std::complex<double> horizontalTerm =
                        receiverHorizontal[mode] * propagation;
                    const std::complex<double> verticalTerm =
                        receiverVertical[mode] * propagation;
                    if (parameters.coherent)
                    {
                        pressure += term;
                        horizontal += horizontalTerm;
                        vertical += verticalTerm;
                    }
                    else
                    {
                        pressure += term * term;
                        horizontal += horizontalTerm * horizontalTerm;
                        vertical += verticalTerm * verticalTerm;
                    }
                }
                if (!parameters.coherent)
                {
                    pressure = std::sqrt(pressure);
                    horizontal = std::sqrt(horizontal);
                    vertical = std::sqrt(vertical);
                }
                const double effectiveRange =
                    result.rangesMetres[range] + parameters.rangeOffsets[receiver];
                if (parameters.sourceType == 'R' &&
                    std::abs(effectiveRange) > std::numeric_limits<double>::min())
                {
                    pressure /= std::sqrt(effectiveRange);
                    horizontal /= std::sqrt(effectiveRange);
                    vertical /= std::sqrt(effectiveRange);
                }
                const std::size_t index =
                    (source * result.receiverDepths.size() + receiver) *
                        result.rangesMetres.size() +
                    range;
                result.values[index] = pressure;
                result.horizontalValues[index] = horizontal;
                result.verticalValues[index] = vertical;
            }
        }
    });
    return result;
}

PressureField evaluateAdiabaticField(
    const ModeFileData &modes,
    const FieldParameters &parameters)
{
    validateFieldParameters(parameters);
    std::vector<const ModeProfileData *> profiles;
    profiles.push_back(&modes);
    for (const ModeProfileData &profile : modes.additionalProfiles)
    {
        profiles.push_back(&profile);
    }
    if (profiles.size() < 2 || parameters.profileRangesKm.size() != profiles.size() ||
        parameters.profileRangesKm.front() != 0.0)
    {
        throw std::invalid_argument("adiabatic field profile inputs are inconsistent");
    }
    for (std::size_t profile = 0; profile < profiles.size(); ++profile)
    {
        if (profiles[profile]->modes.empty() ||
            profiles[profile]->wavenumbers.size() != profiles[profile]->modes.size() ||
            profiles[profile]->depths.size() < 2 ||
            (profile > 0 && parameters.profileRangesKm[profile] <=
                                  parameters.profileRangesKm[profile - 1]))
        {
            throw std::invalid_argument("adiabatic MOD profiles are inconsistent");
        }
    }

    PressureField result;
    result.title = parameters.title.empty() ? modes.title : parameters.title;
    result.frequency = modes.frequency;
    result.sourceDepths = parameters.sourceDepths;
    result.receiverDepths = parameters.receiverDepths;
    result.rangesMetres = parameters.rangesMetres;
    result.values.resize(result.sourceDepths.size() * result.receiverDepths.size() *
                         result.rangesMetres.size());
    result.horizontalValues.resize(result.values.size());
    result.verticalValues.resize(result.values.size());
    result.workerThreadsUsed = static_cast<int>(std::min(
        result.sourceDepths.size(),
        static_cast<std::size_t>(std::max(1, parameters.threadCount))));

    constexpr float piSingle = 3.1415926f;
    const std::complex<float> imaginarySingle(0.0f, 1.0f);
    const std::complex<double> factor = roundedComplex(
        imaginarySingle * std::sqrt(2.0f * piSingle) *
        std::exp(imaginarySingle * (piSingle / 4.0f)));
    const std::complex<double> imaginary(0.0, 1.0);

    std::vector<std::vector<std::complex<double>>> receiverModes(profiles.size());
    std::vector<std::vector<std::complex<double>>> receiverDerivatives(profiles.size());
    for (std::size_t profile = 0; profile < profiles.size(); ++profile)
    {
        const std::size_t modesInProfile = profiles[profile]->modes.size();
        receiverModes[profile].resize(result.receiverDepths.size() * modesInProfile);
        receiverDerivatives[profile].resize(result.receiverDepths.size() * modesInProfile);
        for (std::size_t receiver = 0; receiver < result.receiverDepths.size(); ++receiver)
        {
            for (std::size_t mode = 0; mode < modesInProfile; ++mode)
            {
                const std::size_t index = receiver * modesInProfile + mode;
                receiverModes[profile][index] = roundedComplex(
                    interpolateModeAtDepth(*profiles[profile], mode,
                                           result.receiverDepths[receiver],
                                           modes.frequency));
                receiverDerivatives[profile][index] = roundedComplex(
                    interpolateModeDerivativeAtDepth(*profiles[profile], mode,
                                                     result.receiverDepths[receiver],
                                                     modes.frequency));
            }
        }
    }

    parallelFor(result.sourceDepths.size(),
                static_cast<std::size_t>(std::max(1, parameters.threadCount)),
                [&](std::size_t source) {
        const double sourceSoundSpeed = source < parameters.sourceSoundSpeeds.size()
            ? parameters.sourceSoundSpeeds[source] : 1500.0;
        const double omega = 2.0 * pi * modes.frequency;
        std::size_t modeCount = std::min(
            profiles[0]->modes.size(),
            static_cast<std::size_t>(std::max(0, parameters.modeLimit)));
        modeCount = std::min(modeCount, profiles[1]->modes.size());
        std::vector<std::complex<double>> constants(modeCount);
        for (std::size_t mode = 0; mode < modeCount; ++mode)
        {
            std::complex<double> sourceMode = roundedComplex(
                interpolateModeAtDepth(*profiles[0], mode,
                                       result.sourceDepths[source], modes.frequency));
            if (parameters.beamPattern && source == 0)
            {
                const std::complex<double> k =
                    roundedComplex(profiles[0]->wavenumbers[mode]);
                const double omega = 2.0 * pi * modes.frequency;
                const double kz2 = std::max(
                    0.0, (omega * omega / (1500.0 * 1500.0) - k * k).real());
                const double angle = std::atan(std::sqrt(kz2) / k.real()) *
                                     180.0 / pi;
                sourceMode *= interpolateBeam(parameters, angle);
            }
            constants[mode] = factor * sourceMode;
        }
        std::vector<std::complex<double>> sumK(modeCount, 0.0);
        std::vector<std::complex<double>> sumKInv(modeCount, 0.0);
        std::size_t profile = 0;

        for (std::size_t rangeIndex = 0; rangeIndex < result.rangesMetres.size(); ++rangeIndex)
        {
            const double receiverRange = result.rangesMetres[rangeIndex];
            while (profile + 1 < profiles.size() &&
                   receiverRange > 1000.0 * parameters.profileRangesKm[profile + 1])
            {
                const double left = rangeIndex > 0
                    ? std::max(result.rangesMetres[rangeIndex - 1],
                               1000.0 * parameters.profileRangesKm[profile])
                    : 1000.0 * parameters.profileRangesKm[profile];
                const double right = 1000.0 * parameters.profileRangesKm[profile + 1];
                const double midpoint = 0.5 * (right + left);
                const double weight = (midpoint / 1000.0 -
                                       parameters.profileRangesKm[profile]) /
                                      (parameters.profileRangesKm[profile + 1] -
                                       parameters.profileRangesKm[profile]);
                for (std::size_t mode = 0; mode < modeCount; ++mode)
                {
                    const std::complex<double> kMid = roundedComplex(
                        profiles[profile]->wavenumbers[mode] + weight *
                        (profiles[profile + 1]->wavenumbers[mode] -
                         profiles[profile]->wavenumbers[mode]));
                    sumK[mode] += kMid * (right - left);
                    sumKInv[mode] += (right - left) / kMid;
                }
                ++profile;
                if (profile + 1 < profiles.size())
                {
                    modeCount = std::min(modeCount, profiles[profile + 1]->modes.size());
                }
            }

            const double segmentLeft = 1000.0 * parameters.profileRangesKm[profile];
            const double segmentRight = profile + 1 < profiles.size()
                ? 1000.0 * parameters.profileRangesKm[profile + 1]
                : 1.0e23;
            const double left = rangeIndex > 0
                ? std::max(result.rangesMetres[rangeIndex - 1], segmentLeft)
                : segmentLeft;
            const double midpoint = 0.5 * (receiverRange + left);
            const double weight = (receiverRange - segmentLeft) /
                                  (segmentRight - segmentLeft);
            const double midpointWeight = (midpoint - segmentLeft) /
                                          (segmentRight - segmentLeft);
            const std::size_t rightProfile = std::min(profile + 1, profiles.size() - 1);
            std::vector<std::complex<double>> hank(modeCount);
            std::vector<std::complex<double>> kAtRange(modeCount);
            for (std::size_t mode = 0; mode < modeCount; ++mode)
            {
                const std::complex<double> kMid = roundedComplex(
                    profiles[profile]->wavenumbers[mode] + midpointWeight *
                    (profiles[rightProfile]->wavenumbers[mode] -
                     profiles[profile]->wavenumbers[mode]));
                kAtRange[mode] = roundedComplex(
                    profiles[profile]->wavenumbers[mode] + weight *
                    (profiles[rightProfile]->wavenumbers[mode] -
                     profiles[profile]->wavenumbers[mode]));
                sumK[mode] += kMid * (receiverRange - left);
                sumKInv[mode] += (receiverRange - left) / kMid;
                const std::complex<double> propagation = parameters.coherent
                    ? roundedComplex(std::exp(-imaginary * sumK[mode]))
                    : std::complex<double>(
                          std::exp((-imaginary * sumK[mode]).real()), 0.0);
                hank[mode] = constants[mode] * propagation;
                if (parameters.sourceType == 'R')
                {
                    hank[mode] = receiverRange == 0.0 ? std::complex<double>(0.0, 0.0)
                                                      : hank[mode] /
                                                            std::sqrt(kAtRange[mode] * receiverRange);
                }
                else if (parameters.sourceType == 'X')
                {
                    hank[mode] /= kAtRange[mode];
                }
                else if (parameters.sourceType == 'T')
                {
                    hank[mode] /= std::sqrt(kAtRange[mode] * sumKInv[mode]);
                }
                else
                {
                    hank[mode] /= std::sqrt(kAtRange[mode]);
                }
            }

            for (std::size_t receiver = 0; receiver < result.receiverDepths.size(); ++receiver)
            {
                std::complex<double> pressure = 0.0;
                std::complex<double> horizontal = 0.0;
                std::complex<double> vertical = 0.0;
                double energy = 0.0;
                double horizontalEnergy = 0.0;
                double verticalEnergy = 0.0;
                for (std::size_t mode = 0; mode < modeCount; ++mode)
                {
                    const std::complex<double> leftMode = receiverModes[profile][
                        receiver * profiles[profile]->modes.size() + mode];
                    const std::complex<double> rightMode = receiverModes[rightProfile][
                        receiver * profiles[rightProfile]->modes.size() + mode];
                    const std::complex<double> phi = roundedComplex(
                        leftMode + weight * (rightMode - leftMode));
                    const std::complex<double> term = phi * hank[mode];
                    const std::complex<double> horizontalTerm =
                        sourceSoundSpeed * kAtRange[mode] * term / omega;
                    const std::complex<double> leftDerivative =
                        receiverDerivatives[profile][
                            receiver * profiles[profile]->modes.size() + mode];
                    const std::complex<double> rightDerivative =
                        receiverDerivatives[rightProfile][
                            receiver * profiles[rightProfile]->modes.size() + mode];
                    const std::complex<double> derivative = roundedComplex(
                        leftDerivative + weight * (rightDerivative - leftDerivative));
                    const std::complex<double> verticalTerm =
                        sourceSoundSpeed * derivative * hank[mode] /
                        (omega * imaginary);
                    if (parameters.coherent)
                    {
                        pressure += term;
                        horizontal += horizontalTerm;
                        vertical += verticalTerm;
                    }
                    else
                    {
                        energy += std::norm(term);
                        horizontalEnergy += std::norm(horizontalTerm);
                        verticalEnergy += std::norm(verticalTerm);
                    }
                }
                if (!parameters.coherent)
                {
                    pressure = {std::sqrt(energy), 0.0};
                    horizontal = {std::sqrt(horizontalEnergy), 0.0};
                    vertical = {std::sqrt(verticalEnergy), 0.0};
                }
                const std::size_t index =
                    (source * result.receiverDepths.size() + receiver) *
                        result.rangesMetres.size() + rangeIndex;
                result.values[index] = roundedComplex(pressure);
                result.horizontalValues[index] = roundedComplex(horizontal);
                result.verticalValues[index] = roundedComplex(vertical);
            }
        }
    });
    return result;
}

PressureField evaluateField(
    const ModeFileData &modes,
    const FieldParameters &parameters)
{
    if (parameters.profileRangesKm.size() == 1)
    {
        return evaluateRangeIndependentField(modes, parameters);
    }
    if (parameters.propagationType == 'A')
    {
        return evaluateAdiabaticField(modes, parameters);
    }
    if (parameters.propagationType == 'C')
    {
        return evaluateCoupledField(modes, parameters);
    }
    throw std::invalid_argument("unknown range-dependent propagation type");
}

void writeShadeFile(const PressureField &field,
                    const std::filesystem::path &path)
{
    writeShadeFile(field, path, ShadeDataType::Pressure,
                   Grid_Mode::MODE_R_Rectangular);
}

void writeShadeFile(const PressureField &field,
                    const std::filesystem::path &path,
                    ShadeDataType type,
                    Grid_Mode gridType)
{
    const int nfreq = 1;
    const int ntheta = 1;
    const int nsx = 1;
    const int nsy = 1;
    const int nsz = static_cast<int>(field.sourceDepths.size());
    const int nrz = static_cast<int>(field.receiverDepths.size());
    const int nrr = static_cast<int>(field.rangesMetres.size());
    const std::vector<std::complex<double>> *fieldValues = nullptr;
    switch (type)
    {
    case ShadeDataType::Pressure: fieldValues = &field.values; break;
    case ShadeDataType::VerticalVelocity: fieldValues = &field.verticalValues; break;
    case ShadeDataType::HorizontalVelocity: fieldValues = &field.horizontalValues; break;
    }
    if (nsz < 1 || nrz < 1 || nrr < 1 || !fieldValues ||
        fieldValues->size() != static_cast<std::size_t>(nsz * nrz * nrr))
    {
        throw std::invalid_argument("shade field dimensions are inconsistent");
    }
    const int recordWords = std::max({41, 2 * nfreq, 2 * ntheta, 2 * nsx,
                                      2 * nsy, nsz, nrz, 2 * nrr});
    const std::size_t recordBytes = static_cast<std::size_t>(4 * recordWords);
    std::ofstream stream(path, std::ios::binary | std::ios::trunc);
    if (!stream)
    {
        throw std::runtime_error("unable to create SHD file: " + path.string());
    }
    std::vector<char> record(recordBytes, 0);
    putValue(record, 0, recordWords);
    std::string title = field.title;
    title.resize(80, ' ');
    std::memcpy(record.data() + 4, title.data(), 80);
    writeRecord(stream, 0, record);

    record.assign(recordBytes, 0);
    const char *plotType = gridType == Grid_Mode::MODE_R_Rectangular
                               ? "rectilin  "
                               : "irregular ";
    std::memcpy(record.data(), plotType, 10);
    writeRecord(stream, 1, record);
    record.assign(recordBytes, 0);
    putValue(record, 0, nfreq);
    putValue(record, 4, ntheta);
    putValue(record, 8, nsx);
    putValue(record, 12, nsy);
    putValue(record, 16, nsz);
    putValue(record, 20, nrz);
    putValue(record, 24, nrr);
    putValue(record, 28, field.frequency);
    const double attenuation = 0.0;
    putValue(record, 36, attenuation);
    writeRecord(stream, 2, record);

    const auto writeDoubles = [&](std::size_t index, const std::vector<double> &values) {
        std::vector<char> valuesRecord(recordBytes, 0);
        for (std::size_t value = 0; value < values.size(); ++value)
        {
            putValue(valuesRecord, 8 * value, values[value]);
        }
        writeRecord(stream, index, valuesRecord);
    };
    const auto writeFloats = [&](std::size_t index, const std::vector<double> &values) {
        std::vector<char> valuesRecord(recordBytes, 0);
        for (std::size_t value = 0; value < values.size(); ++value)
        {
            putValue(valuesRecord, 4 * value, static_cast<float>(values[value]));
        }
        writeRecord(stream, index, valuesRecord);
    };
    writeDoubles(3, {field.frequency});
    writeDoubles(4, {0.0});
    writeDoubles(5, {0.0});
    writeDoubles(6, {0.0});
    writeFloats(7, field.sourceDepths);
    writeFloats(8, field.receiverDepths);
    writeDoubles(9, field.rangesMetres);

    std::size_t pressureRecord = 10;
    for (int source = 0; source < nsz; ++source)
    {
        for (int receiver = 0; receiver < nrz; ++receiver)
        {
            record.assign(recordBytes, 0);
            for (int range = 0; range < nrr; ++range)
            {
                const std::complex<double> value = (*fieldValues)[
                    static_cast<std::size_t>((source * nrz + receiver) * nrr + range)];
                putValue(record, static_cast<std::size_t>(8 * range),
                         static_cast<float>(value.real()));
                putValue(record, static_cast<std::size_t>(8 * range + 4),
                         static_cast<float>(value.imag()));
            }
            writeRecord(stream, pressureRecord++, record);
        }
    }
}
}
