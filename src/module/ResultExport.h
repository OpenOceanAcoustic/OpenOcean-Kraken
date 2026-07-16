#ifndef OPEN_OCEAN_KRAKENC_RESULT_EXPORT_H
#define OPEN_OCEAN_KRAKENC_RESULT_EXPORT_H

#include "OpenOceanKrakencParams.h"
#include "algorithm/FieldSolver.h"

#include <filesystem>

namespace OpenOceanKrakenc
{
std::filesystem::path resultPath(const std::string &root, const char *extension);
void exportModeResult(const OOKC_parameters &params,
                      const OOKC_output &output,
                      const std::filesystem::path &path);
void exportShadeResult(const OOKC_parameters &params,
                       const OOKC_output &output,
                       const std::filesystem::path &path,
                       ShadeDataType type);
}

#endif
