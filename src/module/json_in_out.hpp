#ifndef OPEN_OCEAN_KRAKENC_JSON_IN_OUT_HPP
#define OPEN_OCEAN_KRAKENC_JSON_IN_OUT_HPP

#include "OpenOceanKrakencParams.h"
#include "nlohmann/json.hpp"

#include <string>

namespace OpenOceanKrakenc
{
using OpenOcean_json = nlohmann::ordered_json;

void to_json(OpenOcean_json &out, const SSP_Mode &value);
void from_json(const OpenOcean_json &in, SSP_Mode &value);
void to_json(OpenOcean_json &out, const Media_Mode &value);
void from_json(const OpenOcean_json &in, Media_Mode &value);
void to_json(OpenOcean_json &out, const AttenuationUnit &value);
void from_json(const OpenOcean_json &in, AttenuationUnit &value);
void to_json(OpenOcean_json &out, const OceanAbsorptionModel &value);
void from_json(const OpenOcean_json &in, OceanAbsorptionModel &value);
void to_json(OpenOcean_json &out, const BC_Mode &value);
void from_json(const OpenOcean_json &in, BC_Mode &value);
void to_json(OpenOcean_json &out, const Source_Mode &value);
void from_json(const OpenOcean_json &in, Source_Mode &value);
void to_json(OpenOcean_json &out, const Grid_Mode &value);
void from_json(const OpenOcean_json &in, Grid_Mode &value);
void to_json(OpenOcean_json &out, const Run_Mode &value);
void from_json(const OpenOcean_json &in, Run_Mode &value);
void to_json(OpenOcean_json &out, const CoherenceType &value);
void from_json(const OpenOcean_json &in, CoherenceType &value);
void to_json(OpenOcean_json &out, const ModeType &value);
void from_json(const OpenOcean_json &in, ModeType &value);

void to_json(OpenOcean_json &out, const Atten_Mode &value);
void from_json(const OpenOcean_json &in, Atten_Mode &value);
void to_json(OpenOcean_json &out, const HSInfo &value);
void from_json(const OpenOcean_json &in, HSInfo &value);
void to_json(OpenOcean_json &out, const Position &value);
void from_json(const OpenOcean_json &in, Position &value);
void to_json(OpenOcean_json &out, const ReflectionCoef &value);
void from_json(const OpenOcean_json &in, ReflectionCoef &value);
void to_json(OpenOcean_json &out, const InternalReflectionCoefInfo &value);
void from_json(const OpenOcean_json &in, InternalReflectionCoefInfo &value);
void to_json(OpenOcean_json &out, const ReflectionCoefInfo &value);
void from_json(const OpenOcean_json &in, ReflectionCoefInfo &value);
void to_json(OpenOcean_json &out, const SrcBmPat &value);
void from_json(const OpenOcean_json &in, SrcBmPat &value);
void to_json(OpenOcean_json &out, const BdryType &value);
void from_json(const OpenOcean_json &in, BdryType &value);
void to_json(OpenOcean_json &out, const FreqInfo &value);
void from_json(const OpenOcean_json &in, FreqInfo &value);
void to_json(OpenOcean_json &out, const MeshParams &value);
void from_json(const OpenOcean_json &in, MeshParams &value);

namespace ssp
{
void to_json(OpenOcean_json &out, const SSPLayer &value);
void from_json(const OpenOcean_json &in, SSPLayer &value);
void to_json(OpenOcean_json &out, const Range_Independent_Area &value);
void from_json(const OpenOcean_json &in, Range_Independent_Area &value);
void to_json(OpenOcean_json &out, const FlattenedData &value);
void from_json(const OpenOcean_json &in, FlattenedData &value);
void to_json(OpenOcean_json &out, const SSPStructure &value);
void from_json(const OpenOcean_json &in, SSPStructure &value);
}

void to_json(OpenOcean_json &out, const OOKC_parameters &params);
void from_json(const OpenOcean_json &in, OOKC_parameters &params);

bool read_json_file(const std::string &path, OOKC_parameters &params);
bool write_json_file(const std::string &path, const OOKC_parameters &params);
std::string parameters_to_json_string(const OOKC_parameters &params);
}

#endif
