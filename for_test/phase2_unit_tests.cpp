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

    const RootResult limited = complexSecant({10.0, 10.0}, 1.0e-15, 1, polynomial);
    require(!limited.converged && limited.failure == RootFailure::IterationLimit,
            "iteration limit must be reported structurally");

    const std::filesystem::path workspace = OPENOCEANKRAKENC_WORKSPACE_DIR;
    const AcousticCase munk = readAcousticEnv(workspace / "test" / "MunkKleaky.env");
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
