#include "algorithm/ComplexNumerics.h"
#include "algorithm/ComplexRootFinder.h"
#include "algorithm/AcousticCase.h"
#include "algorithm/ComplexMatrixBuilder.h"
#include "algorithm/ComplexDispersion.h"
#include "algorithm/KrakencSolver.h"

#include <cmath>
#include <complex>
#include <filesystem>
#include <fstream>
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

bool near(std::complex<double> actual, std::complex<double> expected, double tolerance)
{
    return std::abs(actual - expected) <= tolerance;
}

AcousticCase closedWaveguideCase()
{
    AcousticCase result;
    result.title = "closed-waveguide-test";
    result.frequency = 50.0;
    result.referenceFrequency = 50.0;
    result.cLow = 1500.0;
    result.cHigh = 5000.0;
    result.rMaxKm = 0.0;
    result.attenuationUnit = 'W';
    result.interpolation = AcousticInterpolation::CLinear;
    result.top.type = AcousticBoundaryType::Vacuum;
    result.bottom.type = AcousticBoundaryType::Rigid;
    AcousticLayer layer;
    layer.baseMesh = 100;
    layer.topDepth = 0.0;
    layer.bottomDepth = 100.0;
    layer.samples = {
        {0.0, 1500.0, 0.0, 1.0, 0.0, 0.0},
        {100.0, 1500.0, 0.0, 1.0, 0.0, 0.0},
    };
    result.layers.push_back(layer);
    return result;
}
}

int main()
{
    require(near(pekerisRoot({-4.0, 0.0}), {0.0, 2.0}, 1.0e-14),
            "Pekeris negative-real branch mismatch");
    require(near(pekerisRoot({3.0, 4.0}), {2.0, 1.0}, 1.0e-14),
            "Pekeris ordinary branch mismatch");

    const auto lossy = complexSoundSpeed(1600.0, 0.8, 50.0, 'W');
    const double expectedImaginary = 0.8 * 1600.0 / (8.6858896 * 2.0 * pi);
    require(std::abs(lossy.real() - 1600.0) <= 1.0e-13,
            "complex sound-speed real part mismatch");
    require(std::abs(lossy.imag() - expectedImaginary) <= 1.0e-12,
            "W attenuation conversion mismatch");

    AttenuationContext powerLaw;
    powerLaw.unit = 'm';
    powerLaw.referenceFrequency = 100.0;
    powerLaw.beta = 2.0;
    powerLaw.transitionFrequency = 200.0;
    for (const double frequency : {50.0, 200.0, 400.0})
    {
        powerLaw.frequency = frequency;
        const double scale = frequency < powerLaw.transitionFrequency
                                 ? std::pow(frequency / powerLaw.referenceFrequency,
                                            powerLaw.beta)
                                 : (frequency / powerLaw.referenceFrequency) *
                                       std::pow(powerLaw.transitionFrequency /
                                                    powerLaw.referenceFrequency,
                                                powerLaw.beta - 1.0);
        const double nepers = 0.8 / 8.6858896 * scale;
        const double expected = nepers * 1600.0 * 1600.0 /
                                (2.0 * pi * frequency);
        require(near(complexSoundSpeed(75.0, 1600.0, 0.8, powerLaw),
                     {1600.0, expected}, 1.0e-12),
                "lowercase-m power-law conversion mismatch");
    }
    AttenuationContext metres = powerLaw;
    metres.unit = 'M';
    metres.frequency = 50.0;
    const double metresExpected = (0.8 / 8.6858896) * 1600.0 * 1600.0 /
                                  (2.0 * pi * metres.frequency);
    require(near(complexSoundSpeed(75.0, 1600.0, 0.8, metres),
                 {1600.0, metresExpected}, 1.0e-12),
            "uppercase M must not use lowercase-m power-law parameters");

    require(std::abs(thorpNepersPerMetre(10000.0) -
                     1.3698423460369063e-4) <= 1.0e-16,
            "Thorp 10 kHz reference mismatch");
    VolumeAbsorptionParameters volume;
    volume.temperatureCelsius = 10.0;
    volume.salinityPsu = 35.0;
    volume.ph = 8.0;
    volume.meanDepthMetres = 1000.0;
    require(std::abs(francoisGarrisonNepersPerMetre(10000.0, volume) -
                     9.732982808053584e-5) <= 1.0e-16,
            "Francois-Garrison Table-IV reference mismatch");
    volume.biologicalLayers.push_back({100.0, 200.0, 1000.0, 5.0, 2.0});
    const double biologicalExpected =
        2.0 / ((1.0 - 1000.0 * 1000.0 / (2000.0 * 2000.0)) *
                   (1.0 - 1000.0 * 1000.0 / (2000.0 * 2000.0)) +
               1.0 / 25.0) /
        8685.8896;
    require(std::abs(biologicalNepersPerMetre(150.0, 2000.0, volume) -
                     biologicalExpected) <= 1.0e-16,
            "biological in-layer absorption mismatch");
    require(biologicalNepersPerMetre(99.0, 2000.0, volume) == 0.0 &&
                biologicalNepersPerMetre(201.0, 2000.0, volume) == 0.0,
            "biological absorption leaked outside its depth layer");

    const ScaledFunction polynomial = [](std::complex<double> z) {
        return ScaledComplex{z * z + std::complex<double>(1.0, 0.0), 0};
    };
    const RootResult root = complexSecant({0.1, 0.9}, 1.0e-13, 80, polynomial);
    require(root.converged, "complex secant did not converge");
    require(root.failure == RootFailure::None, "converged root reported failure");
    require(near(root.root, {0.0, 1.0}, 1.0e-10), "complex secant root mismatch");
    require(root.log10Residual < -10.0, "complex secant residual is too large");

    const RootResult invalid = complexSecant({1.0, 0.0}, 0.0, 10, polynomial);
    require(!invalid.converged && invalid.failure == RootFailure::InvalidTolerance,
            "non-positive tolerance must be rejected structurally");
    require(std::isinf(invalid.relativeCorrection),
            "failed root must not report a zero relative correction");

    const RootResult limited = complexSecant({10.0, 10.0}, 1.0e-15, 1, polynomial);
    require(!limited.converged && limited.failure == RootFailure::IterationLimit,
            "iteration limit must be reported structurally");

    const std::filesystem::path workspace = OPENOCEANKRAKENC_WORKSPACE_DIR;
    const AcousticCase munk = readAcousticEnv(
        workspace / "test" / "toolbox_env" / "MunkKleaky.env");
    require(munk.title == "Munk profile, leaky modes", "Munk title parse mismatch");
    require(std::abs(munk.frequency - 50.0) <= 1.0e-14, "Munk frequency parse mismatch");
    require(munk.layers.size() == 1, "Munk acoustic layer count mismatch");
    require(munk.layers.front().baseMesh == 5000, "Munk base mesh parse mismatch");
    require(munk.layers.front().samples.size() == 27, "Munk SSP sample count mismatch");
    require(munk.interpolation == AcousticInterpolation::N2Linear,
            "Munk interpolation type mismatch");
    require(munk.enableRootRestarts,
            "Munk root-restart option mismatch");
    require(munk.top.type == AcousticBoundaryType::Vacuum,
            "Munk top boundary mismatch");
    require(munk.bottom.type == AcousticBoundaryType::HalfSpace,
            "Munk bottom boundary mismatch");
    require(std::abs(munk.bottom.cp - 1600.0) <= 1.0e-14,
            "Munk bottom sound speed mismatch");
    require(std::abs(munk.bottom.alphaP - 0.8) <= 1.0e-14,
            "Munk bottom attenuation mismatch");

    const std::vector<AcousticCase> solve3Profiles = readAcousticEnvironments(
        workspace / "test" / "toolbox_env" / "solve3_mode_gain.env");
    require(solve3Profiles.size() == 2,
            "solve3 continuation ENV profile count mismatch");
    require(solve3Profiles[0].layers.size() == 1 &&
                solve3Profiles[1].layers.size() == 1 &&
                solve3Profiles[0].layers[0].samples.size() == 2 &&
                solve3Profiles[1].layers[0].samples.size() == 2,
            "solve3 continuation ENV topology mismatch");
    for (std::size_t profile = 0; profile < solve3Profiles.size(); ++profile)
    {
        const AcousticLayer &layer = solve3Profiles[profile].layers.front();
        require(layer.samples.front().cs == 0.0 &&
                    layer.samples.front().rho == 1.0 &&
                    layer.samples.front().alphaP == 0.0 &&
                    layer.samples.front().alphaS == 0.0,
                "solve3 first-row material inheritance mismatch");
        require(layer.samples.back().cs == 0.0 &&
                    layer.samples.back().rho == 1.0 &&
                    layer.samples.back().alphaP == 0.0 &&
                    layer.samples.back().alphaS == 0.0,
                "solve3 within-profile material inheritance mismatch");
    }
    require(solve3Profiles[0].layers.front().samples.front().cp == 1500.0 &&
                solve3Profiles[1].layers.front().samples.front().cp == 1485.0,
            "solve3 profile-specific sound speed was not preserved");
    bool inheritanceMismatchRejected = false;
    try
    {
        static_cast<void>(readAcousticEnvironments(
            std::filesystem::path(OPENOCEANKRAKENC_SOURCE_DIR) /
            "for_test" / "fixtures" / "solve3_inheritance_mismatch.env"));
    }
    catch (const std::runtime_error &error)
    {
        const std::string message = error.what();
        inheritanceMismatchRejected =
            message.find("profile 2") != std::string::npos &&
            message.find("medium 1") != std::string::npos &&
            message.find("depth 0") != std::string::npos;
    }
    require(inheritanceMismatchRejected,
            "mismatched profile inheritance lacks profile/medium/depth context");

    const std::filesystem::path topHalfSpaceEnv =
        std::filesystem::temp_directory_path() / "openocean_krakenc_top_halfspace.env";
    {
        std::ofstream stream(topHalfSpaceEnv);
        stream << "'top acoustic half-space test'\n"
                  "50.0\n"
                  "1\n"
                  "'NAW  '\n"
                  "0.0 1600.0 0.0 1.8 0.8 /\n"
                  "100 0.0 100.0\n"
                  "0.0 1500.0 0.0 1.0 0.0 0.0\n"
                  "100.0 1500.0 /\n"
                  "'R' 0.0\n"
                  "1500.0 5000.0\n"
                  "0.0\n";
    }
    const AcousticCase topHalfSpace = readAcousticEnv(topHalfSpaceEnv);
    std::filesystem::remove(topHalfSpaceEnv);
    require(topHalfSpace.top.type == AcousticBoundaryType::HalfSpace &&
                std::abs(topHalfSpace.top.cp - 1600.0) <= 1.0e-14 &&
                std::abs(topHalfSpace.top.alphaP - 0.8) <= 1.0e-14,
            "top acoustic half-space parse mismatch");

    const auto parseTemporaryEnv = [](
        const std::string &name, const std::string &contents) {
        const std::filesystem::path path =
            std::filesystem::temp_directory_path() / name;
        {
            std::ofstream stream(path);
            stream << contents;
        }
        AcousticCase result = readAcousticEnv(path);
        std::filesystem::remove(path);
        return result;
    };
    const AcousticCase lowercaseM = parseTemporaryEnv(
        "openocean_krakenc_lowercase_m.env",
        "'lowercase m and Thorp'\n"
        "50.0\n1\n'CVmT  '\n"
        "10 0.02 100.0 2.0 80.0\n"
        "0.0 1500.0 0.0 1.0 0.8 0.0\n"
        "100.0 1500.0 /\n"
        "'R' 0.03 1.5 90.0\n"
        "1400.0 2000.0\n0.0\n");
    require(lowercaseM.attenuationUnit == 'm' &&
                lowercaseM.absorptionModel == OceanAbsorptionModel::Thorpe &&
                lowercaseM.layers.front().roughnessRms == 0.02 &&
                lowercaseM.layers.front().attenuationPower == 2.0 &&
                lowercaseM.layers.front().transitionFrequency == 80.0 &&
                lowercaseM.bottom.roughnessRms == 0.03 &&
                lowercaseM.bottom.attenuationPower == 1.5 &&
                lowercaseM.bottom.transitionFrequency == 90.0,
            "lowercase-m/Thorp ENV context parse mismatch");
    const AcousticMatrix lowercaseMatrix = buildAcousticMatrix(lowercaseM, 1);
    const std::complex<double> expectedLowercaseCp = complexSoundSpeed(
        0.0, 1500.0, 0.8,
        attenuationContext(lowercaseM, 2.0, 80.0));
    require(near(lowercaseMatrix.cp.front(), expectedLowercaseCp, 1.0e-12),
            "lowercase-m/Thorp ENV context was not applied to the mesh");

    const AcousticCase francoisGarrison = parseTemporaryEnv(
        "openocean_krakenc_francois_garrison.env",
        "'Francois-Garrison'\n"
        "10000.0\n1\n'CVWF  '\n"
        "10.0 35.0 8.0 1000.0\n"
        "10 0.0 100.0\n"
        "0.0 1500.0 0.0 1.0 0.0 0.0\n"
        "100.0 1500.0 /\n"
        "'R' 0.0\n1400.0 2000.0\n0.0\n");
    require(francoisGarrison.absorptionModel ==
                OceanAbsorptionModel::FrancGarr &&
                francoisGarrison.volumeAbsorption.temperatureCelsius == 10.0 &&
                francoisGarrison.volumeAbsorption.meanDepthMetres == 1000.0,
            "Francois-Garrison ENV parameters were not parsed");

    const AcousticCase biological = parseTemporaryEnv(
        "openocean_krakenc_biological.env",
        "'biological'\n"
        "2000.0\n1\n'CVWB  '\n"
        "1\n100.0 200.0 1000.0 5.0 2.0\n"
        "20 0.0 300.0\n"
        "0.0 1500.0 0.0 1.0 0.0 0.0\n"
        "150.0 1500.0 /\n300.0 1500.0 /\n"
        "'R' 0.0\n1400.0 2000.0\n0.0\n");
    require(biological.absorptionModel == OceanAbsorptionModel::Biological &&
                biological.volumeAbsorption.biologicalLayers.size() == 1 &&
                biological.volumeAbsorption.biologicalLayers.front().qualityFactor == 5.0,
            "biological ENV parameters were not parsed");

    const AcousticMatrix matrix = buildAcousticMatrix(munk, 1);
    require(matrix.depth.size() == 5001, "Munk mesh point count mismatch");
    require(std::abs(matrix.spacing.front() - 1.0) <= 1.0e-14,
            "Munk mesh spacing mismatch");
    require(std::abs(matrix.depth.front()) <= 1.0e-14 &&
                std::abs(matrix.depth.back() - 5000.0) <= 1.0e-14,
            "Munk mesh depth endpoints mismatch");
    const double expectedMidpoint =
        1.0 / std::sqrt(0.5 / (1548.52 * 1548.52) + 0.5 / (1530.29 * 1530.29));
    require(std::abs(matrix.cp[100].real() - expectedMidpoint) <= 1.0e-10,
            "N2 midpoint interpolation mismatch");
    const std::complex<double> expectedB1 =
        -2.0 + matrix.spacing.front() * matrix.spacing.front() * matrix.omega2 /
                   (matrix.cp[100] * matrix.cp[100]);
    require(near(matrix.b1[100], expectedB1, 1.0e-14), "complex B1 mismatch");

    AcousticCase elastic = munk;
    elastic.layers.front().samples.front().cs = 100.0;
    bool elasticRejected = false;
    try
    {
        static_cast<void>(buildAcousticMatrix(elastic, 1));
    }
    catch (const std::invalid_argument &)
    {
        elasticRejected = true;
    }
    require(elasticRejected, "phase-2 matrix builder must reject elastic media");

    const std::complex<double> trialEigenvalue =
        std::pow(matrix.omega / std::complex<double>(1550.0, 0.0), 2);
    const Impedance vacuum = acousticBoundaryImpedance(
        munk.top, true, matrix.omega2, trialEigenvalue,
        munk.attenuationUnit, munk.frequency);
    require(near(vacuum.f, {1.0, 0.0}, 1.0e-14) &&
                near(vacuum.g, {0.0, 0.0}, 1.0e-14),
            "vacuum impedance mismatch");

    AcousticBoundary rigidBoundary;
    rigidBoundary.type = AcousticBoundaryType::Rigid;
    const Impedance rigidBottom = acousticBoundaryImpedance(
        rigidBoundary, false, matrix.omega2, trialEigenvalue, 'W', 50.0);
    const Impedance rigidTop = acousticBoundaryImpedance(
        rigidBoundary, true, matrix.omega2, trialEigenvalue, 'W', 50.0);
    require(near(rigidBottom.f, {0.0, 0.0}, 1.0e-14) &&
                near(rigidBottom.g, {1.0, 0.0}, 1.0e-14),
            "bottom rigid impedance mismatch");
    require(near(rigidTop.g, {-1.0, 0.0}, 1.0e-14),
            "top boundary impedance sign mismatch");

    const Impedance halfSpace = acousticBoundaryImpedance(
        munk.bottom, false, matrix.omega2, trialEigenvalue,
        munk.attenuationUnit, munk.frequency);
    const std::complex<double> bottomCp = complexSoundSpeed(
        munk.bottom.cp, munk.bottom.alphaP, munk.frequency, munk.attenuationUnit);
    const std::complex<double> expectedGamma = pekerisRoot(
        trialEigenvalue - matrix.omega2 / (bottomCp * bottomCp));
    require(near(halfSpace.f, expectedGamma, 1.0e-14),
            "acoustic half-space gamma mismatch");
    require(near(halfSpace.g, {munk.bottom.rho, 0.0}, 1.0e-14),
            "acoustic half-space density impedance mismatch");

    const ComplexDispersion dispersion(munk, matrix);
    const ScaledComplex rawDispersion = dispersion.evaluate(trialEigenvalue, {});
    require(std::isfinite(rawDispersion.value.real()) &&
                std::isfinite(rawDispersion.value.imag()),
            "acoustic recurrence must return a finite scaled value");
    const std::complex<double> previousRoot = trialEigenvalue - 1.0e-3;
    const ScaledComplex deflated = dispersion.evaluate(trialEigenvalue, {previousRoot});
    const std::complex<double> rawActual =
        rawDispersion.value * std::pow(10.0, rawDispersion.power10);
    const std::complex<double> deflatedActual =
        deflated.value * std::pow(10.0, deflated.power10);
    require(near(deflatedActual, rawActual / (trialEigenvalue - previousRoot),
                 1.0e-9 * std::max(1.0, std::abs(deflatedActual))),
            "dispersion root deflation mismatch");

    const auto restart1 = deterministicRestartPoint(1);
    const auto restart2 = deterministicRestartPoint(2);
    require(std::abs(restart1.first - 0.5) <= 1.0e-14 &&
                std::abs(restart1.second - 1.0 / 3.0) <= 1.0e-14,
            "first deterministic restart point mismatch");
    require(std::abs(restart2.first - 0.25) <= 1.0e-14 &&
                std::abs(restart2.second - 2.0 / 3.0) <= 1.0e-14,
            "second deterministic restart point mismatch");

    const std::complex<double> exactExtrapolation(2.0, 3.0);
    const std::complex<double> errorCoefficient(1.0, -2.0);
    const std::vector<int> multipliers{1, 2, 4};
    std::vector<std::complex<double>> approximations;
    for (int multiplier : multipliers)
    {
        approximations.push_back(
            exactExtrapolation + errorCoefficient / static_cast<double>(multiplier * multiplier));
    }
    require(near(richardsonExtrapolate(multipliers, approximations),
                 exactExtrapolation, 1.0e-13),
            "Richardson extrapolation mismatch");

    const AcousticCase closed = closedWaveguideCase();
    const AcousticSolveResult firstSolve = solveAcousticModes(closed);
    const AcousticSolveResult secondSolve = solveAcousticModes(closed);
    require(firstSolve.modes.size() == 6, "closed waveguide mode count mismatch");
    require(secondSolve.modes.size() == firstSolve.modes.size(),
            "deterministic rerun mode count mismatch");
    const double k0 = 2.0 * pi * closed.frequency / 1500.0;
    for (std::size_t mode = 0; mode < firstSolve.modes.size(); ++mode)
    {
        const double kz = (static_cast<double>(mode) + 0.5) * pi / 100.0;
        const double expectedK = std::sqrt(k0 * k0 - kz * kz);
        require(std::abs(firstSolve.modes[mode].wavenumber.real() - expectedK) <= 5.0e-4,
                "closed waveguide analytic wavenumber mismatch");
        require(near(firstSolve.modes[mode].wavenumber,
                     secondSolve.modes[mode].wavenumber, 1.0e-14),
                "acoustic solver rerun is not deterministic");
        if (mode > 0)
        {
            require(std::abs(firstSolve.modes[mode - 1].eigenvalue -
                             firstSolve.modes[mode].eigenvalue) > 1.0e-10,
                    "duplicate root entered accepted modes");
        }
    }

    std::cout << "PHASE2_UNIT_TESTS_OK\n";
}
