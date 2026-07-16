#include "algorithm/InverseIterationComplex.h"
#include "algorithm/AcousticCase.h"
#include "algorithm/ComplexMatrixBuilder.h"
#include "algorithm/ComplexModeSolver.h"
#include "algorithm/ComplexModeNormalization.h"
#include "algorithm/KrakencSolver.h"
#include "algorithm/ElasticCompound.h"
#include "algorithm/ReflectionBoundary.h"
#include "algorithm/ModeFileWriter.h"

#include <Eigen/Dense>

#include <cmath>
#include <complex>
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

AcousticCase closedWaveguideCase()
{
    AcousticCase result;
    result.title = "closed-waveguide-mode-test";
    result.frequency = 50.0;
    result.referenceFrequency = 50.0;
    result.cLow = 1500.0;
    result.cHigh = 5000.0;
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
    const double eigenvalue = 2.0 - std::sqrt(2.0);
    Eigen::VectorXcd diagonal = Eigen::VectorXcd::Constant(3, 2.0 - eigenvalue);
    Eigen::VectorXcd offDiagonal = Eigen::VectorXcd::Constant(2, -1.0);

    const InverseIterationResult realResult =
        inverseIterationComplex(diagonal, offDiagonal, 30);
    require(realResult.converged, "known tridiagonal eigenpair did not converge");
    require(realResult.failure == InverseIterationFailure::None,
            "converged inverse iteration reported failure");
    require(realResult.residual <= 1.0e-10, "inverse iteration residual is too large");
    require(std::abs(std::abs(realResult.vector[0]) - 0.5) <= 1.0e-8 &&
                std::abs(std::abs(realResult.vector[1]) - std::sqrt(0.5)) <= 1.0e-8 &&
                std::abs(std::abs(realResult.vector[2]) - 0.5) <= 1.0e-8,
            "inverse iteration eigenvector shape mismatch");

    const std::complex<double> complexScale(1.0, 0.25);
    const InverseIterationResult complexResult = inverseIterationComplex(
        diagonal * complexScale, offDiagonal * complexScale, 30);
    require(complexResult.converged && complexResult.residual <= 1.0e-10,
            "complex-scaled tridiagonal eigenpair did not converge");
    require(std::abs(std::abs(complexResult.vector[1]) - std::sqrt(0.5)) <= 1.0e-8,
            "complex scaling changed eigenvector shape");

    const InverseIterationResult invalid = inverseIterationComplex(
        Eigen::VectorXcd::Zero(2), Eigen::VectorXcd::Zero(2), 10);
    require(!invalid.converged &&
                invalid.failure == InverseIterationFailure::InvalidDimensions,
            "invalid tridiagonal dimensions must be rejected structurally");

    const InverseIterationResult limited = inverseIterationComplex(
        diagonal, offDiagonal, 0);
    require(!limited.converged &&
                limited.failure == InverseIterationFailure::IterationLimit,
            "inverse iteration limit must be reported structurally");

    const std::filesystem::path workspace = OPENOCEANKRAKENC_WORKSPACE_DIR;
    const AcousticCase elastic = readAcousticEnv(
        workspace / "test" / "elastic_fd_two_layer.env");
    require(elastic.layers.size() == 2 &&
                elastic.layers[1].samples.front().cs == 1300.0,
            "elastic layer was not parsed from ENV");
    require(elastic.bottom.cs == 1600.0,
            "elastic half-space was not parsed from ENV");
    require(!elastic.enableRootRestarts,
            "elastic ENV unexpectedly enabled root restarts");
    const AcousticMatrix elasticMatrix = buildAcousticMatrix(elastic, 1);
    require(elasticMatrix.layerElastic.size() == 2 &&
                !elasticMatrix.layerElastic[0] && elasticMatrix.layerElastic[1],
            "elastic material classification mismatch");
    require(elasticMatrix.cs[elasticMatrix.offsets[1]].real() > 0.0 &&
                std::isfinite(elasticMatrix.b2[elasticMatrix.offsets[1]].real()),
            "elastic compound-matrix coefficients are unavailable");
    const std::complex<double> elasticTrial = std::pow(
        elasticMatrix.omega / std::complex<double>(1700.0, 0.0), 2);
    CompoundState initialState{{1.0, 0.75, -0.25, 0.5, -0.125}};
    int scalePower = 0;
    CompoundState propagated = initialState;
    propagateElasticDown(elasticMatrix, 1, elasticTrial, propagated, scalePower);
    require(compoundStateFinite(propagated),
            "downward elastic compound propagation is non-finite");
    CompoundState repeated = initialState;
    int repeatedPower = 0;
    propagateElasticDown(elasticMatrix, 1, elasticTrial, repeated, repeatedPower);
    double repeatError = 0.0;
    for (std::size_t index = 0; index < initialState.size(); ++index)
    {
        repeatError = std::max(repeatError,
                               std::abs(propagated[index] - repeated[index]));
    }
    require(repeatError == 0.0 && scalePower == repeatedPower,
            "elastic compound propagation is not deterministic");
    const Impedance effectiveBottom = caseBoundaryImpedance(
        elastic, elasticMatrix, false, elasticTrial);
    require(std::isfinite(effectiveBottom.f.real()) &&
                std::isfinite(effectiveBottom.f.imag()) &&
                std::isfinite(effectiveBottom.g.real()) &&
                std::isfinite(effectiveBottom.g.imag()),
            "elastic stack boundary impedance is non-finite");
    const AcousticSolveResult elasticRoots = solveAcousticModes(elastic);
    require(!elasticRoots.modes.empty(), "elastic case root set is empty");
    const std::complex<double> elasticBaseEigenvalue =
        elasticRoots.meshWavenumberSets.front().front() *
        elasticRoots.meshWavenumberSets.front().front();
    const ComplexModeResult elasticMode = solveAcousticMode(
        elastic, elasticMatrix, elasticBaseEigenvalue);
    require(elasticMode.converged && elasticMode.mode.size() == 401,
            "elastic case acoustic pressure mode extraction failed");
    const NormalizedComplexMode normalizedElastic = normalizeAcousticMode(
        elastic, elasticMatrix, elasticBaseEigenvalue,
        elasticMode.turningPoint, elasticMode.mode);
    require(normalizedElastic.mode.allFinite() &&
                std::isfinite(normalizedElastic.groupVelocity),
            "elastic case acoustic pressure normalization is non-finite");

    const AcousticCase brc = readAcousticEnv(
        workspace / "test" / "neggradC_brc.env");
    require(brc.bottom.type == AcousticBoundaryType::ReflectionCoefficient &&
                brc.bottom.reflectionSamples.size() == 91,
            "BRC boundary table was not parsed");
    const AcousticMatrix brcMatrix = buildAcousticMatrix(brc, 1);
    const std::complex<double> brcTrial = std::pow(
        brcMatrix.omega / std::complex<double>(1600.0, 0.0), 2);
    const Impedance brcImpedance = reflectionBoundaryImpedance(
        brc.bottom, false, brcMatrix.omega2, brcTrial,
        brcMatrix.cp.back());
    require(std::abs(brcImpedance.f) <= 1.0e-14 &&
                std::abs(brcImpedance.g - std::complex<double>(1.0, 0.0)) <= 1.0e-14,
            "unit BRC reflection must reduce to a rigid boundary");

    AcousticBoundary ircBoundary;
    ircBoundary.type = AcousticBoundaryType::InternalReflection;
    ircBoundary.internalSamples = {
        {1.0, {1.0, 0.0}, {2.0, 0.0}, 0},
        {2.0, {4.0, 0.0}, {5.0, 0.0}, 0},
        {3.0, {9.0, 0.0}, {10.0, 0.0}, 0},
    };
    const Impedance ircImpedance = reflectionBoundaryImpedance(
        ircBoundary, false, 1.0, {2.5, 0.25}, {1500.0, 0.0});
    require(std::abs(ircImpedance.f - std::complex<double>(6.1875, 1.25)) <= 1.0e-12 &&
                std::abs(ircImpedance.g - std::complex<double>(7.1875, 1.25)) <= 1.0e-12,
            "IRC three-point complex polynomial interpolation mismatch");

    const AcousticCase closed = closedWaveguideCase();
    const AcousticMatrix closedMatrix = buildAcousticMatrix(closed, 1);
    const AcousticSolveResult closedRoots = solveAcousticModes(closed);
    require(!closedRoots.modes.empty(), "closed waveguide roots are unavailable");
    const ComplexModeResult closedMode = solveAcousticMode(
        closed, closedMatrix, closedRoots.modes.front().eigenvalue);
    require(closedMode.converged, "closed waveguide complex mode did not converge");
    require(closedMode.depth.size() == 101 && closedMode.mode.size() == 101,
            "closed waveguide mode grid size mismatch");
    require(closedMode.residual <= 1.0e-8,
            "closed waveguide mode residual is too large");

    Eigen::VectorXd analytic(closedMode.depth.size());
    for (Eigen::Index index = 0; index < analytic.size(); ++index)
    {
        analytic[index] = std::sin(0.5 * std::acos(-1.0) *
                                   closedMode.depth[index] / 100.0);
    }
    analytic.normalize();
    const double shapeCorrelation =
        std::abs(analytic.cast<std::complex<double>>().dot(closedMode.mode));
    require(shapeCorrelation >= 0.999,
            "closed waveguide mode shape does not match the analytic solution");

    const NormalizedComplexMode normalizedClosed = normalizeAcousticMode(
        closed, closedMatrix, closedRoots.modes.front().eigenvalue,
        closedMode.turningPoint, closedMode.mode);
    require(normalizedClosed.mode.allFinite(),
            "closed waveguide normalized mode is non-finite");
    require(std::abs(normalizedClosed.complexNorm -
                     std::complex<double>(1.0, 0.0)) <= 1.0e-10,
            "closed waveguide complex norm is not unity");
    const double closedK = std::sqrt(closedRoots.modes.front().eigenvalue.real());
    const double expectedGroupVelocity = 1500.0 * 1500.0 * closedK /
                                         closedMatrix.omega;
    require(std::abs(normalizedClosed.groupVelocity - expectedGroupVelocity) <=
                1.0e-6 * expectedGroupVelocity,
            "closed waveguide group velocity mismatch");
    require(normalizedClosed.mode[normalizedClosed.turningPoint].real() >= 0.0,
            "normalized mode turning-point phase is not deterministic");

    const AcousticCase munk = readAcousticEnv(workspace / "test" / "MunkKleaky.env");
    require(munk.sourceDepths.size() == 2 && munk.receiverDepths.size() == 1001,
            "Munk source/receiver depth lists were not parsed");
    const AcousticMatrix munkMatrix = buildAcousticMatrix(munk, 1);
    const AcousticSolveResult munkRoots = solveAcousticModes(munk);
    require(munkRoots.modes.size() == 329, "Munk root count changed before mode extraction");
    for (const ModeRoot &root : munkRoots.modes)
    {
        const ComplexModeResult mode =
            solveAcousticMode(munk, munkMatrix, root.eigenvalue);
        require(mode.converged, "a Munk complex mode did not converge");
        require(mode.mode.size() == 5001 && mode.mode.allFinite(),
                "a Munk complex mode is missing or non-finite");
        require(mode.mode.norm() > 0.0 && mode.residual <= 1.0e-7,
                "a Munk complex mode has an invalid residual");
    }

    AcousticCase writable = closed;
    writable.sourceDepths = {25.0};
    writable.receiverDepths = {0.0, 50.0, 100.0};
    const std::filesystem::path modePath =
        std::filesystem::temp_directory_path() / "openocean_krakenc_phase3.mod";
    const ModeFileWriteResult writeResult = writeModeFile(writable, modePath);
    require(writeResult.modeCount == 6 && writeResult.tabulatedDepthCount == 4 &&
                std::filesystem::file_size(modePath) > 0,
            "C++ MOD writer did not produce a complete mode set");

    const std::filesystem::path serialModePath =
        std::filesystem::temp_directory_path() / "openocean_krakenc_phase3_serial.mod";
    const std::filesystem::path parallelModePath =
        std::filesystem::temp_directory_path() / "openocean_krakenc_phase3_parallel.mod";
    const std::vector<AcousticCase> profiles{writable, writable};
    static_cast<void>(writeModeFile(profiles, serialModePath, 1));
    static_cast<void>(writeModeFile(profiles, parallelModePath, 2));
    std::ifstream serialStream(serialModePath, std::ios::binary);
    std::ifstream parallelStream(parallelModePath, std::ios::binary);
    const std::vector<char> serialBytes{
        std::istreambuf_iterator<char>(serialStream), std::istreambuf_iterator<char>()};
    const std::vector<char> parallelBytes{
        std::istreambuf_iterator<char>(parallelStream), std::istreambuf_iterator<char>()};
    require(!serialBytes.empty() && serialBytes == parallelBytes,
            "single-thread and multi-thread MOD files differ");

    AcousticCase allElastic = writable;
    for (AcousticSample &sample : allElastic.layers.front().samples)
    {
        sample.cs = 700.0;
    }
    bool rejectedAllElastic = false;
    try
    {
        static_cast<void>(writeModeFile(allElastic, modePath));
    }
    catch (const std::invalid_argument &)
    {
        rejectedAllElastic = true;
    }
    require(rejectedAllElastic,
            "all-elastic MOD input must fail before matrix or record indexing");
    AcousticCase emptySamples = writable;
    emptySamples.layers.front().samples.clear();
    bool rejectedEmptySamples = false;
    try
    {
        static_cast<void>(writeModeFile(emptySamples, modePath));
    }
    catch (const std::invalid_argument &)
    {
        rejectedEmptySamples = true;
    }
    require(rejectedEmptySamples,
            "empty programmatic layer samples must be rejected before indexing");
    serialStream.close();
    parallelStream.close();
    std::filesystem::remove(modePath);
    std::filesystem::remove(serialModePath);
    std::filesystem::remove(parallelModePath);

    std::cout << "PHASE3_UNIT_TESTS_OK\n";
}
