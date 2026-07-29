#include "algorithm/AcousticCase.h"
#include "algorithm/ComplexMatrixBuilder.h"
#include "algorithm/ComplexModeNormalization.h"
#include "algorithm/ComplexModeSolver.h"
#include "algorithm/ComplexNumerics.h"
#include "algorithm/ComplexRootFinder.h"
#include "algorithm/KrakencSolver.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

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

bool near(std::complex<double> left,
          std::complex<double> right,
          double absoluteTolerance,
          double relativeTolerance)
{
    return std::abs(left - right) <=
           absoluteTolerance +
                relativeTolerance * std::max(std::abs(left), std::abs(right));
}

using TestFunction = void (*)();

struct TestCase
{
    std::string id;
    TestFunction function = nullptr;
};

std::vector<TestCase> &testCases()
{
    static std::vector<TestCase> value;
    return value;
}

class Registration
{
public:
    Registration(const char *id, TestFunction function)
    {
        testCases().push_back({id, function});
    }
};

void testAlg001()
{
    const RootResult constant = complexSecant(
        {2.0, -1.0}, 1.0e-8, 20,
        [](std::complex<double>) {
            return ScaledComplex{{1.0, 0.0}, 0};
        });
    require(!constant.converged, "constant function falsely converged");
    require(constant.failure == RootFailure::DegenerateSecant,
            "constant function must report DegenerateSecant");
    require(constant.log10Residual == 0.0,
            "constant residual must remain log10(1)");

    const RootResult scaledQuadratic = complexSecant(
        {1.5, 0.0}, 1.0e-12, 80,
        [](std::complex<double> z) {
            return ScaledComplex{z * z - 2.0, 20};
        });
    require(scaledQuadratic.converged,
            "scaled quadratic root did not converge");
    require(std::abs(scaledQuadratic.root -
                     std::complex<double>{std::sqrt(2.0), 0.0}) <= 1.0e-12,
            "scaled quadratic root mismatch");
    require(scaledQuadratic.log10Residual > std::log10(1.0e-12),
            "scaled quadratic did not distinguish raw residual from backward error");

    int extremeExponentEvaluations = 0;
    const RootResult extremeExponent = complexSecant(
        {2.0e102, 0.0}, 1.0e100, 20,
        [&extremeExponentEvaluations](std::complex<double>) {
            ++extremeExponentEvaluations;
            if (extremeExponentEvaluations == 1)
            {
                return ScaledComplex{
                    {1.0, 0.0}, std::numeric_limits<int>::min()};
            }
            if (extremeExponentEvaluations == 2)
            {
                return ScaledComplex{
                    {1.0, 0.0}, std::numeric_limits<int>::max()};
            }
            return ScaledComplex{{0.0, 0.0}, 0};
        });
    require(extremeExponent.converged &&
                near(extremeExponent.root, {3.0e102, 0.0}, 1.0e88, 1.0e-14),
            "extreme exponent alignment changed the root");

    const std::complex<double> expected{2.0, 3.0};
    const RootResult linear = complexSecant(
        {2.01, 3.0}, 1.0e-10, 20,
        [expected](std::complex<double> z) {
            return ScaledComplex{z - expected, 0};
        });
    require(linear.converged, "analytic linear root did not converge");
    require(linear.failure == RootFailure::None,
            "analytic root has a failure status");
    require(near(linear.root, expected, 1.0e-12, 1.0e-12),
            "analytic linear root mismatch");

    const RootResult nonFinite = complexSecant(
        {0.0, 0.0}, 1.0e-8, 20,
        [](std::complex<double>) {
            return ScaledComplex{{
                std::numeric_limits<double>::quiet_NaN(), 0.0}, 0};
        });
    require(!nonFinite.converged &&
                nonFinite.failure == RootFailure::NonFiniteValue,
            "non-finite function value was not rejected");
}

const Registration registerAlg001{"ALG-001", &testAlg001};

AcousticCase closedWaveguideCase()
{
    AcousticCase input;
    input.title = "root-mode-remediation";
    input.frequency = 50.0;
    input.referenceFrequency = 50.0;
    input.cLow = 1500.0;
    input.cHigh = 5000.0;
    input.rMaxKm = 0.0;
    input.top.type = AcousticBoundaryType::Vacuum;
    input.bottom.type = AcousticBoundaryType::Rigid;
    AcousticLayer layer;
    layer.baseMesh = 100;
    layer.topDepth = 0.0;
    layer.bottomDepth = 100.0;
    layer.samples = {
        {0.0, 1500.0, 0.0, 1.0, 0.0, 0.0},
        {100.0, 1500.0, 0.0, 1.0, 0.0, 0.0}};
    input.layers.push_back(layer);
    return input;
}

void testAlg002()
{
    const RootUniquenessSpec spec{};
    require(sameEigenroot(
                {{1.0 + 5.0e-15, 0.0}, 1.0e-14},
                {{1.0, 0.0}, 1.0e-14}, 1.0e-3, spec),
            "true duplicate not detected");
    require(!sameEigenroot(
                {{1.0 + 1.97392088021789e-9, 0.0}, 1.0e-14},
                {{1.0, 0.0}, 1.0e-14},
                1.97392088021789e-9, spec),
            "separate local modes were merged");

    // RootCandidateGate is the same state machine used by solveMesh. This
    // scripted sequence forces duplicate -> advance -> distinct acceptance.
    RootCandidateGate gate(spec, 16);
    require(gate.submit({{1.0, 0.0}, 1.0e-14}) ==
                RootCandidateAction::Accepted,
            "first scripted root was not accepted");
    require(gate.submit({{1.0 + 5.0e-15, 0.0}, 1.0e-14}) ==
                RootCandidateAction::DuplicateAdvance,
            "forced duplicate did not advance the candidate loop");
    require(gate.accepted().size() == 1 &&
                gate.consecutiveDuplicateAdvances() == 1,
            "forced duplicate corrupted candidate-loop state");
    require(gate.submit({{1.0 - 1.97392088021789e-9, 0.0}, 1.0e-14}) ==
                RootCandidateAction::Accepted,
            "candidate loop did not continue after the forced duplicate");
    require(gate.accepted().size() == 2 &&
                gate.consecutiveDuplicateAdvances() == 0,
            "post-duplicate acceptance did not reset the counter");

    int overflowingNumeratorEvaluations = 0;
    const RootResult overflowingNumerator = complexSecant(
        {1.0, 0.0}, 1.0, 20,
        [&overflowingNumeratorEvaluations](std::complex<double>) {
            ++overflowingNumeratorEvaluations;
            if (overflowingNumeratorEvaluations == 1)
            {
                return ScaledComplex{{0.0, 0.0}, 0};
            }
            if (overflowingNumeratorEvaluations == 2)
            {
                return ScaledComplex{{
                    std::numeric_limits<double>::max() / 2.0, 0.0}, 0};
            }
            return ScaledComplex{{0.0, 0.0}, 0};
        });
    require(overflowingNumerator.converged &&
                near(overflowingNumerator.root, {0.9, 0.0}, 1.0e-12, 1.0e-12),
            "overflowing secant numerator bypassed krakenc limiting");

    AcousticCase input = closedWaveguideCase();
    input.frequency = input.referenceFrequency = 0.1;
    input.layers.front().baseMesh = 200;
    input.layers.front().bottomDepth = 100000.0;
    input.layers.front().samples.back().depth = 100000.0;

    const double k0 = 2.0 * pi * 0.1 / 1500.0;
    const double lambda0 =
        k0 * k0 - std::pow(pi / 200000.0, 2);
    const double lambda1 =
        k0 * k0 - std::pow(3.0 * pi / 200000.0, 2);
    AcousticSolveResult result;
    try
    {
        result = solveAcousticModes(input);
    }
    catch (const std::exception &)
    {
        require(false,
                "deep-water source-aligned convergence was rejected");
    }
    require(result.modes.size() >= 2,
            "resolvable deep-water modes merged");
    require(std::abs(result.modes[0].eigenvalue.real() - lambda0) <=
                5.0e-12,
            "first analytic eigenvalue mismatch");
    require(std::abs(result.modes[1].eigenvalue.real() - lambda1) <=
                5.0e-12,
            "second analytic eigenvalue mismatch");
    require(std::any_of(
                result.rootSearchDiagnostics.begin(),
                result.rootSearchDiagnostics.end(),
                [](const RootSearchDiagnostic &diagnostic) {
                    return diagnostic.accepted &&
                           diagnostic.log10Residual > -14.0;
                }),
            "accepted deep-water root did not expose the raw residual plateau");
}

const Registration registerAlg002{"ALG-002", &testAlg002};

void testAlg006()
{
    std::vector<ModeRoot> comparisonModes(3);
    comparisonModes[0].eigenvalue = {9.0, 0.0};
    comparisonModes[1].eigenvalue = {6.0, 0.0};
    comparisonModes[2].eigenvalue = {3.0, 8.0};
    require(krakencFinalModeCount(comparisonModes, 2.0, 1.0) == 2,
            "ALG-006 final cutoff did not use real extrapolated eigenvalue");

    std::vector<ModeRoot> boundaryMode(1);
    boundaryMode[0].eigenvalue = {4.0, 8.0};
    require(krakencFinalModeCount(boundaryMode, 2.0, 1.0) == 0,
            "ALG-006 final cutoff did not exclude the strict boundary");

    std::vector<ModeRoot> prefixModes(4);
    prefixModes[0].eigenvalue = {9.0, 0.0};
    prefixModes[1].eigenvalue = {3.0, 0.0};
    prefixModes[2].eigenvalue = {6.0, 0.0};
    prefixModes[3].eigenvalue = {2.0, 0.0};
    require(krakencFinalModeCount(prefixModes, 2.0, 1.0) == 3,
            "ALG-006 final cutoff did not preserve MINLOC prefix semantics");

    std::vector<ModeRoot> scatteredMode(1);
    scatteredMode[0].eigenvalue = {9.0, 0.0};
    scatteredMode[0].scatterPerturbation = {-8.0, 0.0};
    scatteredMode[0].wavenumber = {1.0, 0.0};
    require(krakencFinalModeCount(scatteredMode, 2.0, 1.0) == 1,
            "ALG-006 final cutoff incorrectly used scatter-adjusted wavenumber");
}

const Registration registerAlg006{"ALG-006", &testAlg006};

void testAlg009()
{
    constexpr double scatterAbsoluteTolerance = 1.0e-15;
    constexpr double scatterRelativeTolerance = 1.0e-12;
    constexpr double eigenvalueAbsoluteTolerance = 1.0e-14;
    constexpr double eigenvalueRelativeTolerance = 1.0e-12;
    constexpr double wavenumberAbsoluteTolerance = 1.0e-14;
    constexpr double wavenumberRelativeTolerance = 1.0e-12;

    AcousticCase smooth = closedWaveguideCase();
    smooth.layers.front().baseMesh = 20;
    smooth.rMaxKm = 1.0;
    const AcousticSolveResult smoothResult = solveAcousticModes(smooth);
    require(smoothResult.meshSetsUsed >= 2,
            "ALG-009 fixture did not use multiple meshes");
    require(smoothResult.modes.size() >= 2,
            "ALG-009 fixture did not retain two modes");
    require(!smoothResult.meshModeCounts.empty(),
            "ALG-009 fixture did not record mesh mode counts");
    const int preCutoffCommonTracks = *std::min_element(
        smoothResult.meshModeCounts.begin(),
        smoothResult.meshModeCounts.end());
    require(preCutoffCommonTracks >= 0 &&
                static_cast<std::size_t>(preCutoffCommonTracks) ==
                    smoothResult.modes.size(),
            "ALG-009 fixture pre-cutoff and retained counts differed");

    std::vector<ModeRoot> probes;
    for (const ModeRoot &mode : smoothResult.modes)
    {
        const bool newCoarseIndex = std::none_of(
            probes.begin(), probes.end(),
            [&](const ModeRoot &probe) {
                return probe.coarseModeIndex == mode.coarseModeIndex;
            });
        if (newCoarseIndex)
        {
            probes.push_back(mode);
        }
        if (probes.size() == 2)
        {
            break;
        }
    }
    require(probes.size() == 2 &&
                probes[0].coarseModeIndex != probes[1].coarseModeIndex,
            "ALG-009 fixture did not provide two coarse mode indices");
    require(probes[0].coarseModeIndex != 0 ||
                probes[1].coarseModeIndex != 0,
            "ALG-009 fixture did not provide a nonzero coarse mode index");
    for (const ModeRoot &probe : probes)
    {
        require(probe.hasCoarseProvenance,
                "ALG-009 fixture probe lacked coarse provenance");
        require(std::abs(probe.coarseEigenvalue - probe.eigenvalue) > 1.0e-12,
                "ALG-009 fixture coarse and extrapolated roots were too close");
    }

    const auto finalizeWavenumber = [](std::complex<double> value) {
        std::complex<double> result = std::sqrt(value);
        if (result.imag() > 0.0)
        {
            result = {result.real(), 0.0};
        }
        return result;
    };

    struct RoughCaseOracle
    {
        AcousticCase rough;
        std::vector<std::complex<double>> expectedScatter;
    };

    const auto buildRoughOracle = [&](const AcousticCase &rough) {
        const AcousticMatrix coarse = buildAcousticMatrix(rough, 1);
        for (const ModeRoot &mode : smoothResult.modes)
        {
            const ComplexModeResult rawCoarse = solveAcousticMode(
                rough, coarse, mode.coarseEigenvalue);
            require(rawCoarse.converged,
                    "ALG-009 coarse inverse-iteration fixture did not converge");
        }

        std::vector<std::complex<double>> expectedScatter;
        expectedScatter.reserve(probes.size());
        for (const ModeRoot &probe : probes)
        {
            const ComplexModeResult rawCoarse = solveAcousticMode(
                rough, coarse, probe.coarseEigenvalue);
            require(rawCoarse.converged,
                    "ALG-009 coarse probe inverse iteration did not converge");
            const NormalizedComplexMode expected = normalizeAcousticMode(
                rough, coarse, probe.coarseEigenvalue,
                rawCoarse.turningPoint, rawCoarse.mode);
            require(std::isfinite(expected.scatterPerturbation.real()) &&
                        std::isfinite(expected.scatterPerturbation.imag()),
                    "ALG-009 coarse scatter oracle was not finite");
            require(std::abs(expected.scatterPerturbation) >
                        100.0 * scatterAbsoluteTolerance,
                    "ALG-009 coarse scatter oracle lacked magnitude");
            require(!near(
                        probe.coarseEigenvalue, probe.eigenvalue,
                        eigenvalueAbsoluteTolerance,
                        eigenvalueRelativeTolerance),
                    "ALG-009 fixture roots were equal at fixed tolerance");

            const std::complex<double> correctWavenumber =
                finalizeWavenumber(
                    probe.eigenvalue + expected.scatterPerturbation);
            const std::complex<double> wrongWavenumber =
                finalizeWavenumber(
                    probe.coarseEigenvalue + expected.scatterPerturbation);
            const double wavenumberEnvelope =
                wavenumberAbsoluteTolerance +
                wavenumberRelativeTolerance *
                    std::max(std::abs(correctWavenumber),
                             std::abs(wrongWavenumber));
            require(!near(
                        correctWavenumber, wrongWavenumber,
                        wavenumberAbsoluteTolerance,
                        wavenumberRelativeTolerance) &&
                        std::abs(correctWavenumber - wrongWavenumber) >
                            100.0 * wavenumberEnvelope,
                    "ALG-009 final-wavenumber fixture lacked discrimination");

            const ComplexModeResult rawMixed = solveAcousticMode(
                rough, coarse, probe.eigenvalue);
            bool mixedDiscriminates = !rawMixed.converged;
            if (rawMixed.converged)
            {
                const NormalizedComplexMode mixed = normalizeAcousticMode(
                    rough, coarse, probe.eigenvalue,
                    rawMixed.turningPoint, rawMixed.mode);
                const double scatterEnvelope =
                    scatterAbsoluteTolerance +
                    scatterRelativeTolerance *
                        std::max(
                            std::abs(expected.scatterPerturbation),
                            std::abs(mixed.scatterPerturbation));
                mixedDiscriminates =
                    !near(
                        mixed.scatterPerturbation,
                        expected.scatterPerturbation,
                        scatterAbsoluteTolerance,
                        scatterRelativeTolerance) &&
                    std::abs(
                        mixed.scatterPerturbation -
                        expected.scatterPerturbation) >
                        100.0 * scatterEnvelope;
            }
            require(mixedDiscriminates,
                    "ALG-009 mixed-root negative control lacked discrimination");
            expectedScatter.push_back(expected.scatterPerturbation);
        }
        return RoughCaseOracle{rough, expectedScatter};
    };

    struct RoughCaseResult
    {
        bool completed = false;
        AcousticSolveResult result;
    };

    const auto solveRoughCase = [](const AcousticCase &rough) {
        RoughCaseResult outcome;
        try
        {
            outcome.result = solveAcousticModes(rough);
            outcome.completed = true;
        }
        catch (const std::runtime_error &error)
        {
            if (std::string(error.what()) !=
                "roughness mode extraction failed before scatter correction")
            {
                throw;
            }
        }
        return outcome;
    };

    const auto verifyRoughCase = [&](const RoughCaseOracle &oracle,
                                     const RoughCaseResult &outcome) {
        for (std::size_t index = 0; index < probes.size(); ++index)
        {
            const ModeRoot &probe = probes[index];
            const auto actual = std::find_if(
                outcome.result.modes.begin(), outcome.result.modes.end(),
                [&](const ModeRoot &mode) {
                    return mode.coarseModeIndex == probe.coarseModeIndex;
                });
            const bool paired =
                outcome.completed &&
                actual != outcome.result.modes.end() &&
                near(
                    actual->scatterPerturbation,
                    oracle.expectedScatter[index],
                    scatterAbsoluteTolerance,
                    scatterRelativeTolerance);
            require(
                paired,
                "ALG-009 coarse matrix/eigenvalue pairing contract was not satisfied");
            require(near(
                        actual->eigenvalue, probe.eigenvalue,
                        eigenvalueAbsoluteTolerance,
                        eigenvalueRelativeTolerance),
                    "ALG-009 scattering replaced extrapolated eigenvalue");
            const std::complex<double> expectedWavenumber =
                finalizeWavenumber(
                    probe.eigenvalue + oracle.expectedScatter[index]);
            require(near(
                        actual->wavenumber, expectedWavenumber,
                        wavenumberAbsoluteTolerance,
                        wavenumberRelativeTolerance),
                    "ALG-009 final wavenumber did not combine extrapolated eigenvalue with coarse scatter");
        }
    };

    AcousticCase bottomRough = smooth;
    bottomRough.bottom.roughnessRms = 0.25;
    AcousticCase layerRough = smooth;
    layerRough.layers.front().roughnessRms = 0.25;

    const RoughCaseOracle bottomOracle = buildRoughOracle(bottomRough);
    const RoughCaseOracle layerOracle = buildRoughOracle(layerRough);
    const RoughCaseResult bottomResult = solveRoughCase(bottomRough);
    const RoughCaseResult layerResult = solveRoughCase(layerRough);
    verifyRoughCase(bottomOracle, bottomResult);
    verifyRoughCase(layerOracle, layerResult);
}

const Registration registerAlg009{"ALG-009", &testAlg009};
}

int main(int argc, char **argv)
{
    if (argc != 2)
    {
        return 2;
    }
    const std::string id = argv[1];
    const auto found = std::find_if(
        testCases().begin(), testCases().end(),
        [&](const TestCase &test) { return test.id == id; });
    if (found == testCases().end() || found->function == nullptr)
    {
        return 2;
    }
    found->function();
    return 0;
}
