#ifndef OPEN_OCEAN_KRAKENC_PARAMETER_ADAPTERS_H
#define OPEN_OCEAN_KRAKENC_PARAMETER_ADAPTERS_H

#include "OpenOceanKrakencParams.h"
#include "algorithm/AcousticCase.h"
#include "algorithm/FieldSolver.h"

#include <vector>

namespace OpenOceanKrakenc
{
std::vector<AcousticCase> toAcousticCases(const OOKC_parameters &params);
FieldParameters toFieldParameters(const OOKC_parameters &params);
void updatePublicParameters(const std::vector<AcousticCase> &cases,
                            OOKC_parameters &params);
void updatePublicParameters(const FieldParameters &field,
                            OOKC_parameters &params);
void validatePublicParameters(const OOKC_parameters &params,
                              Run_Mode requestedRun);
void validatePublicParameterSemantics(const OOKC_parameters &params,
                                      Run_Mode requestedRun);
void validatePublicExecutionCapabilities(const OOKC_parameters &params,
                                         Run_Mode requestedRun);
}

#endif
