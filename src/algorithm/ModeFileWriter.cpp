#include "algorithm/ModeFileWriter.h"

#include "algorithm/ComplexMatrixBuilder.h"
#include "algorithm/ComplexModeNormalization.h"
#include "algorithm/ComplexModeSolver.h"
#include "algorithm/ComplexNumerics.h"
#include "algorithm/KrakencSolver.h"

#include <algorithm>
#include <array>
#include <atomic>
#include <complex>
#include <cstring>
#include <fstream>
#include <future>
#include <stdexcept>
#include <string>

namespace OpenOceanKrakenc
{
namespace
{
template <typename Value>
void put(std::vector<char> &record, std::size_t offset, const Value &value)
{
    if (offset + sizeof(Value) > record.size())
    {
        throw std::runtime_error("MOD record overflow");
    }
    std::memcpy(record.data() + offset, &value, sizeof(Value));
}

void putBytes(std::vector<char> &record, std::size_t offset,
              const char *data, std::size_t count)
{
    if (offset + count > record.size())
    {
        throw std::runtime_error("MOD record overflow");
    }
    std::memcpy(record.data() + offset, data, count);
}

void writeRecord(std::ofstream &stream, std::size_t zeroBasedRecord,
                 const std::vector<char> &record)
{
    static_cast<void>(zeroBasedRecord);
    stream.write(record.data(), static_cast<std::streamsize>(record.size()));
    if (!stream)
    {
        throw std::runtime_error("failed while writing MOD record");
    }
}

char boundaryCode(AcousticBoundaryType type)
{
    switch (type)
    {
    case AcousticBoundaryType::Vacuum:
        return 'V';
    case AcousticBoundaryType::Rigid:
        return 'R';
    case AcousticBoundaryType::HalfSpace:
        return 'A';
    case AcousticBoundaryType::ReflectionCoefficient:
        return 'F';
    case AcousticBoundaryType::InternalReflection:
        return 'P';
    }
    return 'V';
}

std::vector<double> mergedDepths(const AcousticCase &input)
{
    std::vector<double> result = input.sourceDepths;
    result.insert(result.end(), input.receiverDepths.begin(), input.receiverDepths.end());
    if (result.empty())
    {
        for (const AcousticLayer &layer : input.layers)
        {
            if (layer.samples.front().cs == 0.0)
            {
                result.push_back(layer.topDepth);
                result.push_back(layer.bottomDepth);
            }
        }
    }
    std::sort(result.begin(), result.end());
    result.erase(std::unique(result.begin(), result.end(), [](double left, double right) {
                     return std::abs(left - right) <= 1.0e-8;
                 }),
                 result.end());
    return result;
}

void putComplexFloat(std::vector<char> &record, std::size_t offset,
                     std::complex<double> value)
{
    struct ComplexFloat
    {
        float real;
        float imaginary;
    };
    const ComplexFloat output{static_cast<float>(value.real()),
                              static_cast<float>(value.imag())};
    put(record, offset, output);
}

template <typename Function>
void parallelFor(std::size_t count, std::size_t requestedThreads,
                 const Function &function)
{
    const std::size_t workerCount = std::min(
        std::max<std::size_t>(1, requestedThreads), count);
    if (workerCount == 1)
    {
        for (std::size_t index = 0; index < count; ++index)
        {
            function(index);
        }
        return;
    }
    std::atomic<std::size_t> next{0};
    std::vector<std::future<void>> workers;
    workers.reserve(workerCount);
    for (std::size_t worker = 0; worker < workerCount; ++worker)
    {
        workers.push_back(std::async(std::launch::async, [&]() {
            for (;;)
            {
                const std::size_t index = next.fetch_add(1);
                if (index >= count)
                {
                    return;
                }
                function(index);
            }
        }));
    }
    for (auto &worker : workers)
    {
        worker.get();
    }
}
}

std::vector<ModeFileWriteResult> writeModeFile(
    const std::vector<AcousticCase> &inputs,
    const std::filesystem::path &outputPath,
    std::size_t threadCount)
{
    if (inputs.empty())
    {
        throw std::invalid_argument("cannot write an empty MOD profile set");
    }
    for (const AcousticCase &input : inputs)
    {
        if (input.layers.empty() || std::any_of(
                input.layers.begin(), input.layers.end(), [](const AcousticLayer &layer) {
                    return layer.samples.empty();
                }))
        {
            throw std::invalid_argument("MOD output requires nonempty layer samples");
        }
        const bool hasAcousticLayer = std::any_of(
            input.layers.begin(), input.layers.end(), [](const AcousticLayer &layer) {
                return !layer.samples.empty() && layer.samples.front().cs == 0.0;
            });
        if (!hasAcousticLayer)
        {
            throw std::invalid_argument("MOD output requires at least one acoustic layer");
        }
    }
    int recordWords = 32;
    for (const AcousticCase &input : inputs)
    {
        const int depthCount = static_cast<int>(mergedDepths(input).size());
        const int acousticMedia = static_cast<int>(std::count_if(
            input.layers.begin(), input.layers.end(), [](const AcousticLayer &layer) {
                return !layer.samples.empty() && layer.samples.front().cs == 0.0;
            }));
        recordWords = std::max({recordWords, 2 * depthCount, 3 * acousticMedia});
    }
    const std::size_t recordBytes = static_cast<std::size_t>(4 * recordWords);
    std::ofstream stream(outputPath, std::ios::binary | std::ios::trunc);
    if (!stream)
    {
        throw std::runtime_error("unable to create MOD file: " + outputPath.string());
    }
    std::vector<ModeFileWriteResult> results;
    std::vector<AcousticSolveResult> solves(inputs.size());
    parallelFor(inputs.size(), threadCount, [&](std::size_t index) {
        solves[index] = solveAcousticModes(inputs[index]);
    });

    std::size_t profileRecord = 0;
    for (std::size_t inputIndex = 0; inputIndex < inputs.size(); ++inputIndex)
    {
    const AcousticCase &input = inputs[inputIndex];
    const AcousticSolveResult &solve = solves[inputIndex];
    if (solve.modes.empty() || solve.meshWavenumberSets.empty() ||
        solve.meshWavenumberSets.front().size() < solve.modes.size())
    {
        throw std::runtime_error("base mesh does not contain every output mode");
    }
    const AcousticMatrix matrix = buildAcousticMatrix(input, 1);
    const std::vector<double> depths = mergedDepths(input);

    std::size_t firstAcoustic = input.layers.size();
    std::size_t lastAcoustic = input.layers.size();
    for (std::size_t medium = 0; medium < matrix.layerElastic.size(); ++medium)
    {
        if (!matrix.layerElastic[medium])
        {
            if (firstAcoustic == input.layers.size())
            {
                firstAcoustic = medium;
            }
            lastAcoustic = medium;
        }
    }
    if (firstAcoustic == input.layers.size())
    {
        throw std::invalid_argument("MOD output requires at least one acoustic layer");
    }
    const int acousticMedia = static_cast<int>(lastAcoustic - firstAcoustic + 1);
    const int modeCount = static_cast<int>(solve.modes.size());
    const int depthCount = static_cast<int>(depths.size());
    std::vector<char> record(recordBytes, 0);
    put(record, 0, recordWords);
    std::string title = "KRAKENC " + input.title;
    title.resize(80, ' ');
    putBytes(record, 4, title.data(), 80);
    const int one = 1;
    put(record, 84, one);
    put(record, 88, acousticMedia);
    put(record, 92, depthCount);
    put(record, 96, depthCount);
    writeRecord(stream, profileRecord, record);

    record.assign(recordBytes, 0);
    std::size_t cursor = 0;
    for (std::size_t medium = firstAcoustic; medium <= lastAcoustic; ++medium)
    {
        put(record, cursor, matrix.counts[medium]);
        cursor += 4;
        static constexpr char material[8] = {'A', 'C', 'O', 'U', 'S', 'T', 'I', 'C'};
        putBytes(record, cursor, material, 8);
        cursor += 8;
    }
    writeRecord(stream, profileRecord + 1, record);

    record.assign(recordBytes, 0);
    cursor = 0;
    for (std::size_t medium = firstAcoustic; medium <= lastAcoustic; ++medium)
    {
        const float topDepth = static_cast<float>(input.layers[medium].topDepth);
        const float rho = static_cast<float>(matrix.rho[matrix.offsets[medium]]);
        put(record, cursor, topDepth);
        put(record, cursor + 4, rho);
        cursor += 8;
    }
    writeRecord(stream, profileRecord + 2, record);

    record.assign(recordBytes, 0);
    put(record, 0, input.frequency);
    writeRecord(stream, profileRecord + 3, record);

    record.assign(recordBytes, 0);
    for (int index = 0; index < depthCount; ++index)
    {
        const float value = static_cast<float>(depths[static_cast<std::size_t>(index)]);
        put(record, static_cast<std::size_t>(4 * index), value);
    }
    writeRecord(stream, profileRecord + 4, record);

    record.assign(recordBytes, 0);
    put(record, 0, modeCount);
    writeRecord(stream, profileRecord + 5, record);

    record.assign(recordBytes, 0);
    cursor = 0;
    const auto putBoundary = [&](const AcousticBoundary &boundary, double depth) {
        put(record, cursor, boundaryCode(boundary.type));
        cursor += 1;
        const std::complex<double> cp = boundary.cp > 0.0
                                            ? complexSoundSpeed(boundary.cp, boundary.alphaP,
                                                                input.frequency, input.attenuationUnit)
                                            : std::complex<double>(0.0, 0.0);
        const std::complex<double> cs = boundary.cs > 0.0
                                            ? complexSoundSpeed(boundary.cs, boundary.alphaS,
                                                                input.frequency, input.attenuationUnit)
                                            : std::complex<double>(0.0, 0.0);
        putComplexFloat(record, cursor, cp);
        cursor += 8;
        putComplexFloat(record, cursor, cs);
        cursor += 8;
        put(record, cursor, static_cast<float>(boundary.rho));
        cursor += 4;
        put(record, cursor, static_cast<float>(depth));
        cursor += 4;
    };
    putBoundary(input.top, input.layers.front().topDepth);
    putBoundary(input.bottom, input.layers.back().bottomDepth);
    writeRecord(stream, profileRecord + 6, record);

    ModeFileWriteResult result;
    result.modeCount = modeCount;
    result.tabulatedDepthCount = depthCount;
    result.recordLengthWords = recordWords;
    struct PreparedMode
    {
        double groupVelocity = 0.0;
        std::vector<char> record;
    };
    std::vector<PreparedMode> preparedModes(static_cast<std::size_t>(modeCount));
    const std::size_t extractionThreads = modeCount >= 4 ? threadCount : 1;
    parallelFor(static_cast<std::size_t>(modeCount), extractionThreads,
                [&](std::size_t preparedIndex) {
        const int modeIndex = static_cast<int>(preparedIndex);
        const std::complex<double> baseK =
            solve.meshWavenumberSets.front()[static_cast<std::size_t>(modeIndex)];
        const std::complex<double> baseEigenvalue = baseK * baseK;
        const ComplexModeResult raw = solveAcousticMode(
            input, matrix, baseEigenvalue);
        if (!raw.converged)
        {
            throw std::runtime_error("inverse iteration failed while writing MOD file");
        }
        const NormalizedComplexMode normalized = normalizeAcousticMode(
            input, matrix, baseEigenvalue, raw.turningPoint, raw.mode);
        PreparedMode &prepared = preparedModes[preparedIndex];
        prepared.groupVelocity = normalized.groupVelocity;
        prepared.record.assign(recordBytes, 0);
        Eigen::Index rightIndex = 1;
        for (int depthIndex = 0; depthIndex < depthCount; ++depthIndex)
        {
            const double target = depths[static_cast<std::size_t>(depthIndex)];
            std::complex<double> value;
            if (target <= raw.depth[0])
            {
                value = normalized.mode[0];
            }
            else if (target >= raw.depth[raw.depth.size() - 1])
            {
                value = normalized.mode[normalized.mode.size() - 1];
            }
            else
            {
                while (rightIndex + 1 < raw.depth.size() &&
                       raw.depth[rightIndex] <= target)
                {
                    ++rightIndex;
                }
                const Eigen::Index leftIndex = rightIndex - 1;
                const double fraction = (target - raw.depth[leftIndex]) /
                                        (raw.depth[rightIndex] - raw.depth[leftIndex]);
                value = (1.0 - fraction) * normalized.mode[leftIndex] +
                        fraction * normalized.mode[rightIndex];
            }
            putComplexFloat(prepared.record,
                            static_cast<std::size_t>(8 * depthIndex), value);
        }
    });
    for (int modeIndex = 0; modeIndex < modeCount; ++modeIndex)
    {
        PreparedMode &prepared = preparedModes[static_cast<std::size_t>(modeIndex)];
        result.groupVelocities.push_back(prepared.groupVelocity);
        writeRecord(stream, profileRecord + static_cast<std::size_t>(7 + modeIndex),
                    prepared.record);
    }

    const int complexPerRecord = recordWords / 2;
    int written = 0;
    int waveRecord = 0;
    while (written < modeCount)
    {
        record.assign(recordBytes, 0);
        const int count = std::min(complexPerRecord, modeCount - written);
        for (int index = 0; index < count; ++index)
        {
            putComplexFloat(record, static_cast<std::size_t>(8 * index),
                            solve.modes[static_cast<std::size_t>(written + index)].wavenumber);
        }
        writeRecord(stream, profileRecord +
                                static_cast<std::size_t>(modeCount + 7 + waveRecord),
                    record);
        written += count;
        ++waveRecord;
    }
    profileRecord += static_cast<std::size_t>(7 + modeCount + waveRecord);
    results.push_back(std::move(result));
    }
    return results;
}

ModeFileWriteResult writeModeFile(
    const AcousticCase &input,
    const std::filesystem::path &outputPath)
{
    std::vector<AcousticCase> inputs{input};
    std::vector<ModeFileWriteResult> results = writeModeFile(inputs, outputPath);
    return std::move(results.front());
}

void writeModeFile(const ModeFileData &modes,
                   const std::filesystem::path &outputPath)
{
    std::vector<const ModeProfileData *> profiles{&modes};
    for (const ModeProfileData &profile : modes.additionalProfiles)
    {
        profiles.push_back(&profile);
    }
    int recordWords = 32;
    for (const ModeProfileData *profile : profiles)
    {
        if (profile->depths.size() < 2 || profile->wavenumbers.empty() ||
            profile->modes.size() != profile->wavenumbers.size() ||
            profile->meshCounts.size() != profile->materials.size() ||
            profile->meshCounts.size() != profile->mediumDepths.size() ||
            profile->meshCounts.size() != profile->mediumDensities.size())
        {
            throw std::invalid_argument("in-memory MOD profile is incomplete");
        }
        for (const auto &mode : profile->modes)
        {
            if (mode.size() != profile->depths.size())
            {
                throw std::invalid_argument("in-memory MOD mode depth count mismatch");
            }
        }
        recordWords = std::max(
            {recordWords,
             static_cast<int>(2 * profile->depths.size()),
             static_cast<int>(3 * profile->meshCounts.size())});
    }
    const std::size_t recordBytes = static_cast<std::size_t>(4 * recordWords);
    std::ofstream stream(outputPath, std::ios::binary | std::ios::trunc);
    if (!stream)
    {
        throw std::runtime_error("unable to create MOD file: " + outputPath.string());
    }
    std::size_t profileRecord = 0;
    for (const ModeProfileData *profile : profiles)
    {
        const int mediumCount = static_cast<int>(profile->meshCounts.size());
        const int depthCount = static_cast<int>(profile->depths.size());
        const int modeCount = static_cast<int>(profile->wavenumbers.size());
        std::vector<char> record(recordBytes, 0);
        put(record, 0, recordWords);
        std::string title = "KRAKENC " + modes.title;
        title.resize(80, ' ');
        putBytes(record, 4, title.data(), 80);
        const int one = 1;
        put(record, 84, one);
        put(record, 88, mediumCount);
        put(record, 92, depthCount);
        put(record, 96, depthCount);
        writeRecord(stream, profileRecord, record);

        record.assign(recordBytes, 0);
        std::size_t cursor = 0;
        for (int medium = 0; medium < mediumCount; ++medium)
        {
            put(record, cursor, profile->meshCounts[static_cast<std::size_t>(medium)]);
            cursor += 4;
            std::string material = profile->materials[static_cast<std::size_t>(medium)];
            material.resize(8, ' ');
            putBytes(record, cursor, material.data(), 8);
            cursor += 8;
        }
        writeRecord(stream, profileRecord + 1, record);

        record.assign(recordBytes, 0);
        cursor = 0;
        for (int medium = 0; medium < mediumCount; ++medium)
        {
            put(record, cursor, static_cast<float>(profile->mediumDepths[static_cast<std::size_t>(medium)]));
            put(record, cursor + 4, static_cast<float>(profile->mediumDensities[static_cast<std::size_t>(medium)]));
            cursor += 8;
        }
        writeRecord(stream, profileRecord + 2, record);

        record.assign(recordBytes, 0);
        put(record, 0, modes.frequency);
        writeRecord(stream, profileRecord + 3, record);

        record.assign(recordBytes, 0);
        for (int depth = 0; depth < depthCount; ++depth)
        {
            put(record, static_cast<std::size_t>(4 * depth),
                static_cast<float>(profile->depths[static_cast<std::size_t>(depth)]));
        }
        writeRecord(stream, profileRecord + 4, record);

        record.assign(recordBytes, 0);
        put(record, 0, modeCount);
        writeRecord(stream, profileRecord + 5, record);

        record.assign(recordBytes, 0);
        cursor = 0;
        const auto putBoundaryData = [&](const ModeBoundaryData &boundary) {
            put(record, cursor, boundary.type);
            cursor += 1;
            putComplexFloat(record, cursor, boundary.cp);
            cursor += 8;
            putComplexFloat(record, cursor, boundary.cs);
            cursor += 8;
            put(record, cursor, static_cast<float>(boundary.rho));
            cursor += 4;
            put(record, cursor, static_cast<float>(boundary.depth));
            cursor += 4;
        };
        putBoundaryData(profile->top);
        putBoundaryData(profile->bottom);
        writeRecord(stream, profileRecord + 6, record);

        for (int mode = 0; mode < modeCount; ++mode)
        {
            record.assign(recordBytes, 0);
            for (int depth = 0; depth < depthCount; ++depth)
            {
                putComplexFloat(record, static_cast<std::size_t>(8 * depth),
                                profile->modes[static_cast<std::size_t>(mode)]
                                              [static_cast<std::size_t>(depth)]);
            }
            writeRecord(stream, profileRecord + static_cast<std::size_t>(7 + mode), record);
        }

        const int complexPerRecord = recordWords / 2;
        int written = 0;
        int waveRecord = 0;
        while (written < modeCount)
        {
            record.assign(recordBytes, 0);
            const int count = std::min(complexPerRecord, modeCount - written);
            for (int index = 0; index < count; ++index)
            {
                putComplexFloat(record, static_cast<std::size_t>(8 * index),
                                profile->wavenumbers[static_cast<std::size_t>(written + index)]);
            }
            writeRecord(stream, profileRecord +
                                    static_cast<std::size_t>(modeCount + 7 + waveRecord),
                        record);
            written += count;
            ++waveRecord;
        }
        profileRecord += static_cast<std::size_t>(7 + modeCount + waveRecord);
    }
}
}
