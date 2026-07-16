#include "algorithm/FieldSolver.h"
#include "algorithm/AcousticCase.h"

#include <cmath>
#include <complex>
#include <filesystem>
#include <iostream>
#include <stdexcept>

using namespace OpenOceanKrakenc;

namespace
{
void require(bool condition, const char *message)
{
    if (!condition)
    {
        throw std::runtime_error(message);
    }
}
}

int main()
{
    const std::filesystem::path workspace = OPENOCEANKRAKENC_WORKSPACE_DIR;
    const FieldParameters munk = readFieldParameters(
        workspace / "test" / "MunkKleaky.flp");
    require(munk.sourceType == 'R' && munk.beamPattern && munk.coherent,
            "Munk FLP options parse mismatch");
    require(munk.rangesMetres.size() == 501 &&
                std::abs(munk.rangesMetres.back() - 100000.0) <= 1.0e-9,
            "Munk receiver range expansion mismatch");
    require(munk.sourceDepths.size() == 2 && munk.receiverDepths.size() == 1001 &&
                munk.rangeOffsets.size() == munk.receiverDepths.size(),
            "Munk FLP depth/offset parse mismatch");
    const FieldParameters omni = readFieldParameters(
        workspace / "OpenOcean-Krakenc" / "for_test" / "fixtures" /
        "omni_small.flp");
    require(omni.sourceType == 'R' && !omni.beamPattern,
            "legal FLP omni flag O was not accepted");

    const std::vector<AcousticCase> stepProfiles = readAcousticEnvironments(
        workspace / "test" / "stepK_rd.env");
    require(stepProfiles.size() == 4,
            "stepK multi-profile ENV count mismatch");
    require(std::abs(soundSpeedAt(stepProfiles.front(), 18.0) - 1476.7) <= 1.0e-9,
            "stepK source sound-speed interpolation mismatch");
    const std::vector<AcousticCase> wedgeProfiles = readAcousticEnvironments(
        workspace / "test" / "wedge.env");
    require(wedgeProfiles.size() == 51,
            "wedge multi-profile ENV count mismatch");

    ModeFileData modeData;
    modeData.title = "single mode analytic field";
    modeData.frequency = 50.0;
    modeData.depths = {0.0, 100.0};
    modeData.top.depth = 0.0;
    modeData.bottom.depth = 100.0;
    modeData.wavenumbers = {{1.0, 0.0}};
    modeData.modes = {{{1.0, 0.0}, {2.0, 0.0}}};

    FieldParameters analytic;
    analytic.sourceType = 'S';
    analytic.coherent = true;
    analytic.modeLimit = 1;
    analytic.profileRangesKm = {0.0};
    analytic.rangesMetres = {10.0};
    analytic.sourceDepths = {50.0};
    analytic.receiverDepths = {100.0};
    analytic.rangeOffsets = {0.0};
    const PressureField pressure = evaluateRangeIndependentField(modeData, analytic);
    require(pressure.values.size() == 1, "analytic field dimensions mismatch");
    const std::complex<double> imaginary(0.0, 1.0);
    const std::complex<double> factor =
        imaginary * std::sqrt(2.0 * std::acos(-1.0)) *
        std::exp(imaginary * std::acos(-1.0) / 4.0);
    const std::complex<double> expected =
        factor * 3.0 * std::exp(-imaginary * 10.0);
    require(std::abs(pressure.values.front() - expected) <= 1.0e-12,
            "single-mode analytic field mismatch");
    require(pressure.horizontalValues.size() == 1 &&
                pressure.verticalValues.size() == 1,
            "analytic velocity field dimensions mismatch");
    const double omega = 2.0 * std::acos(-1.0) * modeData.frequency;
    const std::complex<double> expectedHorizontal = expected * 1500.0 / omega;
    const std::complex<double> expectedVertical =
        expected * 1500.0 * 0.01 / (omega * imaginary * 2.0);
    require(std::abs(pressure.horizontalValues.front() - expectedHorizontal) <= 1.0e-12,
            "single-mode analytic horizontal velocity mismatch");
    require(std::abs(pressure.verticalValues.front() - expectedVertical) <= 1.0e-12,
            "single-mode analytic vertical velocity mismatch");

    ModeFileData halfspaceModes = modeData;
    halfspaceModes.top.type = 'A';
    halfspaceModes.top.cp = {1400.0, 0.0};
    halfspaceModes.top.rho = 1.5;
    FieldParameters boundaryFieldParameters = analytic;
    boundaryFieldParameters.sourceDepths = {0.0};
    boundaryFieldParameters.receiverDepths = {0.0};
    const PressureField boundaryField =
        evaluateRangeIndependentField(halfspaceModes, boundaryFieldParameters);
    FieldParameters upperHalfspaceParameters = boundaryFieldParameters;
    upperHalfspaceParameters.sourceDepths = {-10.0};
    upperHalfspaceParameters.receiverDepths = {-10.0};
    const PressureField upperHalfspaceField =
        evaluateRangeIndependentField(halfspaceModes, upperHalfspaceParameters);
    constexpr float piMode = 3.141592f;
    const double halfspaceK =
        static_cast<double>(2.0f * piMode * static_cast<float>(modeData.frequency)) /
        1400.0;
    const float roundedGamma = static_cast<float>(
        std::sqrt(1.0 - halfspaceK * halfspaceK));
    const std::complex<double> expectedHalfspacePressure =
        boundaryField.values.front() * std::exp(-20.0 * roundedGamma);
    require(std::abs(upperHalfspaceField.values.front() -
                     expectedHalfspacePressure) <= 1.0e-10,
            "range-independent halfspace mode extrapolation mismatch");

    FieldParameters differential = analytic;
    differential.rangesMetres = {9.999, 10.0, 10.001};
    differential.receiverDepths = {49.9, 50.0, 50.1};
    differential.rangeOffsets = {0.0, 0.0, 0.0};
    const PressureField differentialField =
        evaluateRangeIndependentField(modeData, differential);
    const std::size_t centre = 4;
    const std::complex<double> pressureRangeDerivative =
        (differentialField.values[5] - differentialField.values[3]) / 0.002;
    const std::complex<double> pressureDepthDerivative =
        (differentialField.values[7] - differentialField.values[1]) / 0.2;
    require(std::abs(differentialField.horizontalValues[centre] -
                     imaginary * 1500.0 * pressureRangeDerivative / omega) <= 1.0e-5,
            "horizontal velocity does not match the pressure range derivative");
    require(std::abs(differentialField.verticalValues[centre] -
                     1500.0 * pressureDepthDerivative / (omega * imaginary)) <= 1.0e-5,
            "vertical velocity does not match the pressure depth derivative");

    ModeFileData adiabaticModes = modeData;
    adiabaticModes.additionalProfiles.push_back(
        static_cast<const ModeProfileData &>(modeData));
    FieldParameters adiabatic = analytic;
    adiabatic.propagationType = 'A';
    adiabatic.profileRangesKm = {0.0, 1.0};
    const PressureField adiabaticField = evaluateField(adiabaticModes, adiabatic);
    require(adiabaticField.horizontalValues.size() == adiabaticField.values.size() &&
                adiabaticField.verticalValues.size() == adiabaticField.values.size(),
            "adiabatic velocity field dimensions mismatch");
    require(std::abs(adiabaticField.horizontalValues.front() - expectedHorizontal) <= 1.0e-5 &&
                std::abs(adiabaticField.verticalValues.front() - expectedVertical) <= 1.0e-5,
            "adiabatic analytic velocity mismatch");

    FieldParameters shadedAdiabatic = adiabatic;
    shadedAdiabatic.beamPattern = true;
    shadedAdiabatic.beamAnglesDegrees = {-90.0, 90.0};
    shadedAdiabatic.beamAmplitudes = {0.5, 0.5};
    const PressureField shadedAdiabaticField =
        evaluateField(adiabaticModes, shadedAdiabatic);
    require(std::abs(shadedAdiabaticField.values.front() -
                     0.5 * adiabaticField.values.front()) <= 1.0e-5,
            "adiabatic source beam pattern was not applied");

    FieldParameters scaledCoupled = adiabatic;
    scaledCoupled.propagationType = 'C';
    scaledCoupled.sourceType = 'S';
    const PressureField scaledCoupledField =
        evaluateField(adiabaticModes, scaledCoupled);
    require(std::abs(scaledCoupledField.values.front() - expected) <= 1.0e-5,
            "coupled scaled-cylindrical source phase mismatch");
    require(std::abs(scaledCoupledField.horizontalValues.front() -
                     expectedHorizontal) <= 1.0e-5 &&
                std::abs(scaledCoupledField.verticalValues.front() -
                         expectedVertical) <= 1.0e-5,
            "coupled scaled-cylindrical velocity mismatch");

    FieldParameters coupled = adiabatic;
    coupled.propagationType = 'C';
    coupled.sourceType = 'R';
    const PressureField coupledField = evaluateField(adiabaticModes, coupled);
    require(coupledField.horizontalValues.size() == coupledField.values.size() &&
                coupledField.verticalValues.size() == coupledField.values.size(),
            "coupled velocity field dimensions mismatch");
    require(std::isfinite(coupledField.horizontalValues.front().real()) &&
                std::isfinite(coupledField.verticalValues.front().real()),
            "coupled analytic velocity is not finite");

    FieldParameters threaded = analytic;
    threaded.sourceDepths = {50.0, 50.0};
    threaded.threadCount = 2;
    const PressureField threadedField = evaluateField(modeData, threaded);
    threaded.threadCount = 1;
    const PressureField sequentialField = evaluateField(modeData, threaded);
    require(threadedField.workerThreadsUsed == 2,
            "field solver did not use the requested source workers");
    require(threadedField.values == sequentialField.values &&
                threadedField.horizontalValues == sequentialField.horizontalValues &&
                threadedField.verticalValues == sequentialField.verticalValues,
            "single-thread and multi-thread fields differ");

    FieldParameters threadedAdiabatic = adiabatic;
    threadedAdiabatic.sourceDepths = {50.0, 50.0};
    threadedAdiabatic.threadCount = 2;
    const PressureField parallelAdiabatic =
        evaluateField(adiabaticModes, threadedAdiabatic);
    threadedAdiabatic.threadCount = 1;
    const PressureField serialAdiabatic =
        evaluateField(adiabaticModes, threadedAdiabatic);
    require(parallelAdiabatic.workerThreadsUsed == 2 &&
                parallelAdiabatic.values == serialAdiabatic.values &&
                parallelAdiabatic.horizontalValues == serialAdiabatic.horizontalValues &&
                parallelAdiabatic.verticalValues == serialAdiabatic.verticalValues,
            "adiabatic single-thread and multi-thread fields differ");

    FieldParameters threadedCoupled = scaledCoupled;
    threadedCoupled.sourceDepths = {50.0, 50.0};
    threadedCoupled.threadCount = 2;
    const PressureField parallelCoupled =
        evaluateField(adiabaticModes, threadedCoupled);
    threadedCoupled.threadCount = 1;
    const PressureField serialCoupled =
        evaluateField(adiabaticModes, threadedCoupled);
    require(parallelCoupled.workerThreadsUsed == 2 &&
                parallelCoupled.values == serialCoupled.values &&
                parallelCoupled.horizontalValues == serialCoupled.horizontalValues &&
                parallelCoupled.verticalValues == serialCoupled.verticalValues,
            "coupled single-thread and multi-thread fields differ");

    ModeFileData emptyModes = modeData;
    emptyModes.wavenumbers.clear();
    emptyModes.modes.clear();
    bool rejectedEmptyModes = false;
    try
    {
        static_cast<void>(evaluateField(emptyModes, analytic));
    }
    catch (const std::invalid_argument &)
    {
        rejectedEmptyModes = true;
    }
    require(rejectedEmptyModes, "empty MOD mode sets must be rejected");

    FieldParameters invalidField = analytic;
    invalidField.sourceType = 'Q';
    bool rejectedSourceType = false;
    try
    {
        static_cast<void>(evaluateField(modeData, invalidField));
    }
    catch (const std::invalid_argument &)
    {
        rejectedSourceType = true;
    }
    require(rejectedSourceType, "unknown FLP source type must be rejected");

    invalidField = analytic;
    invalidField.sourceType = 'T';
    bool rejectedTranslationalType = false;
    try
    {
        static_cast<void>(evaluateField(modeData, invalidField));
    }
    catch (const std::invalid_argument &)
    {
        rejectedTranslationalType = true;
    }
    require(rejectedTranslationalType,
            "unsupported translational source type must be rejected");

    invalidField = analytic;
    invalidField.rangesMetres = {10.0, 5.0};
    bool rejectedDescendingRanges = false;
    try
    {
        static_cast<void>(evaluateField(modeData, invalidField));
    }
    catch (const std::invalid_argument &)
    {
        rejectedDescendingRanges = true;
    }
    require(rejectedDescendingRanges,
            "descending receiver ranges must be rejected before marching");

    std::cout << "PHASE4_UNIT_TESTS_OK\n";
}
