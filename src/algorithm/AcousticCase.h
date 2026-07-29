#ifndef OPEN_OCEAN_KRAKENC_ACOUSTIC_CASE_H
#define OPEN_OCEAN_KRAKENC_ACOUSTIC_CASE_H

#include "OpenOceanKrakencEnums.h"

#include <filesystem>
#include <complex>
#include <string>
#include <vector>

namespace OpenOceanKrakenc
{
enum class AcousticInterpolation
{
    N2Linear,
    CLinear,
    Pchip,
    CubicSpline
};

enum class AcousticBoundaryType
{
    Vacuum,
    Rigid,
    HalfSpace,
    ReflectionCoefficient,
    InternalReflection
};

struct ReflectionSample
{
    double angleDegrees = 0.0;
    double magnitude = 0.0;
    double phaseRadians = 0.0;
};

struct InternalReflectionSample
{
    double eigenvalue = 0.0;
    std::complex<double> f{};
    std::complex<double> g{};
    int power10 = 0;
};

struct AcousticSample
{
    double depth = 0.0;
    double cp = 0.0;
    double cs = 0.0;
    double rho = 0.0;
    double alphaP = 0.0;
    double alphaS = 0.0;
};

struct AcousticLayer
{
    int baseMesh = 0;
    double topDepth = 0.0;
    double bottomDepth = 0.0;
    double roughnessRms = 0.0;
    double attenuationPower = 0.0;
    double transitionFrequency = 0.0;
    std::vector<AcousticSample> samples;
};

struct AcousticBoundary
{
    AcousticBoundaryType type = AcousticBoundaryType::Vacuum;
    double depth = 0.0;
    double cp = 0.0;
    double cs = 0.0;
    double rho = 0.0;
    double alphaP = 0.0;
    double alphaS = 0.0;
    double roughnessRms = 0.0;
    double attenuationPower = 0.0;
    double transitionFrequency = 0.0;
    std::vector<ReflectionSample> reflectionSamples;
    std::vector<InternalReflectionSample> internalSamples;
};

struct AcousticCase
{
    std::string title;
    double frequency = 0.0;
    double referenceFrequency = 0.0;
    double cLow = 0.0;
    double cHigh = 0.0;
    double rMaxKm = 0.0;
    char attenuationUnit = 'W';
    OceanAbsorptionModel absorptionModel = OceanAbsorptionModel::None;
    VolumeAbsorptionParameters volumeAbsorption;
    bool enableRootRestarts = false;
    AcousticInterpolation interpolation = AcousticInterpolation::CLinear;
    AcousticBoundary top;
    AcousticBoundary bottom;
    std::vector<AcousticLayer> layers;
    std::vector<double> sourceDepths;
    std::vector<double> receiverDepths;
};

AttenuationContext attenuationContext(const AcousticCase &input,
                                      double attenuationPower = 0.0,
                                      double transitionFrequency = 0.0);

AcousticCase readAcousticEnv(const std::filesystem::path &path);
std::vector<AcousticCase> readAcousticEnvironments(
    const std::filesystem::path &path);
double soundSpeedAt(const AcousticCase &input, double depth);
}

#endif
