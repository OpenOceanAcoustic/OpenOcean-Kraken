#ifndef OPEN_OCEAN_KRAKENC_FIELD_SOLVER_H
#define OPEN_OCEAN_KRAKENC_FIELD_SOLVER_H

#include "OpenOceanKrakencEnums.h"

#include <complex>
#include <filesystem>
#include <string>
#include <vector>

namespace OpenOceanKrakenc
{
struct FieldParameters
{
    std::string title;
    char sourceType = 'R';
    char propagationType = ' ';
    bool coherent = true;
    bool beamPattern = false;
    int modeLimit = 20000;
    int threadCount = 1;
    std::vector<double> profileRangesKm;
    std::vector<double> rangesMetres;
    std::vector<double> sourceDepths;
    std::vector<double> sourceSoundSpeeds;
    std::vector<double> receiverDepths;
    std::vector<double> rangeOffsets;
    std::vector<double> beamAnglesDegrees;
    std::vector<double> beamAmplitudes;
};

struct ModeBoundaryData
{
    char type = 'V';
    std::complex<double> cp{};
    std::complex<double> cs{};
    double rho = 0.0;
    double depth = 0.0;
};

struct ModeProfileData
{
    std::vector<double> depths;
    std::vector<std::complex<double>> wavenumbers;
    std::vector<std::vector<std::complex<double>>> modes;
    std::vector<int> meshCounts;
    std::vector<std::string> materials;
    std::vector<double> mediumDepths;
    std::vector<double> mediumDensities;
    ModeBoundaryData top;
    ModeBoundaryData bottom;
};

struct ModeFileData : ModeProfileData
{
    std::string title;
    double frequency = 0.0;
    std::vector<ModeProfileData> additionalProfiles;
};

struct PressureField
{
    std::string title;
    double frequency = 0.0;
    std::vector<double> sourceDepths;
    std::vector<double> receiverDepths;
    std::vector<double> rangesMetres;
    std::vector<std::complex<double>> values;
    std::vector<std::complex<double>> horizontalValues;
    std::vector<std::complex<double>> verticalValues;
    int workerThreadsUsed = 1;
};

enum class ShadeDataType
{
    Pressure = 1,
    VerticalVelocity = 2,
    HorizontalVelocity = 3
};

FieldParameters readFieldParameters(const std::filesystem::path &path);
void validateFieldParameters(const FieldParameters &parameters);
ModeFileData readModeFile(const std::filesystem::path &path);
PressureField evaluateRangeIndependentField(
    const ModeFileData &modes,
    const FieldParameters &parameters);
PressureField evaluateAdiabaticField(
    const ModeFileData &modes,
    const FieldParameters &parameters);
PressureField evaluateCoupledField(
    const ModeFileData &modes,
    const FieldParameters &parameters);
PressureField evaluateField(
    const ModeFileData &modes,
    const FieldParameters &parameters);
void writeShadeFile(const PressureField &field,
                    const std::filesystem::path &path);
void writeShadeFile(const PressureField &field,
                    const std::filesystem::path &path,
                    ShadeDataType type,
                    Grid_Mode gridType);
}

#endif
