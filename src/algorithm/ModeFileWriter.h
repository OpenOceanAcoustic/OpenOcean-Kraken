#ifndef OPEN_OCEAN_KRAKENC_MODE_FILE_WRITER_H
#define OPEN_OCEAN_KRAKENC_MODE_FILE_WRITER_H

#include "algorithm/AcousticCase.h"
#include "algorithm/FieldSolver.h"

#include <filesystem>
#include <vector>

namespace OpenOceanKrakenc
{
struct ModeFileWriteResult
{
    int modeCount = 0;
    int tabulatedDepthCount = 0;
    int recordLengthWords = 0;
    std::vector<double> groupVelocities;
};

ModeFileWriteResult writeModeFile(
    const AcousticCase &input,
    const std::filesystem::path &outputPath);
std::vector<ModeFileWriteResult> writeModeFile(
    const std::vector<AcousticCase> &inputs,
    const std::filesystem::path &outputPath,
    std::size_t threadCount = 1);
void writeModeFile(const ModeFileData &modes,
                   const std::filesystem::path &outputPath);
}

#endif
