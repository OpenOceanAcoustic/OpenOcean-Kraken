#include "OpenOceanKrakencKernelInterface.h"

#include "algorithm/AcousticCase.h"
#include "algorithm/ComplexMatrixBuilder.h"
#include "algorithm/ComplexModeNormalization.h"
#include "algorithm/ComplexModeSolver.h"
#include "algorithm/ComplexNumerics.h"
#include "algorithm/FieldSolver.h"
#include "algorithm/KrakencSolver.h"
#include "algorithm/ModeFileWriter.h"
#include "module/env_in_out.hpp"
#include "module/json_in_out.hpp"
#include "module/ParameterAdapters.h"
#include "module/ResultExport.h"
#include "ThreadPool.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <filesystem>
#include <future>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <thread>
#include <utility>

namespace OpenOceanKrakenc
{
namespace
{
struct SampledMode
{
    std::complex<double> value{};
    std::complex<double> derivative{};
};

SampledMode sampleMode(const Eigen::VectorXd &depths,
                       const Eigen::VectorXcd &mode,
                       double target)
{
    if (depths.size() < 2 || mode.size() != depths.size())
    {
        throw std::invalid_argument("interface mode interpolation dimensions mismatch");
    }
    if (target > depths[depths.size() - 1])
    {
        return {};
    }
    Eigen::Index right = 1;
    if (target >= depths[depths.size() - 1])
    {
        right = depths.size() - 1;
    }
    else if (target > depths[0])
    {
        while (right + 1 < depths.size() && depths[right] < target)
        {
            ++right;
        }
    }
    const Eigen::Index left = right - 1;
    const double interval = depths[right] - depths[left];
    if (!(interval > 0.0))
    {
        throw std::runtime_error("interface mode depths must be strictly increasing");
    }
    const double weight = std::clamp((target - depths[left]) / interval, 0.0, 1.0);
    SampledMode result;
    result.value = (1.0 - weight) * mode[left] + weight * mode[right];
    result.derivative = (mode[right] - mode[left]) / interval;
    return result;
}

EigenParams solveInterfaceProfile(const AcousticCase &input)
{
    const AcousticSolveResult result = solveAcousticModes(input);

    EigenParams eigen;
    const int modalCapacity = static_cast<int>(result.modes.size());
    eigen.resize(modalCapacity, result.meshSetsUsed,
                 static_cast<int>(input.sourceDepths.size()),
                 static_cast<int>(input.receiverDepths.size()));
    eigen.M = modalCapacity;
    for (std::size_t set = 0; set < result.meshWavenumberSets.size(); ++set)
    {
        const std::size_t count = std::min(
            result.modes.size(), result.meshWavenumberSets[set].size());
        for (std::size_t mode = 0; mode < count; ++mode)
        {
            const Eigen::Index index = static_cast<Eigen::Index>(
                set * result.modes.size() + mode);
            const std::complex<double> meshK = result.meshWavenumberSets[set][mode];
            eigen.EVMat[index] = meshK * meshK;
            eigen.Extrap[index] = result.modes[mode].eigenvalue;
        }
    }
    const AcousticMatrix matrix = buildAcousticMatrix(input, 1);
    for (Eigen::Index mode = 0; mode < eigen.k.size(); ++mode)
    {
        const ModeRoot &source = result.modes[static_cast<std::size_t>(mode)];
        eigen.k[mode] = source.wavenumber;
        const std::complex<double> baseK =
            result.meshWavenumberSets.front()[static_cast<std::size_t>(mode)];
        const std::complex<double> baseEigenvalue = baseK * baseK;
        const ComplexModeResult raw = solveAcousticMode(
            input, matrix, baseEigenvalue);
        if (!raw.converged)
        {
            throw std::runtime_error("interface complex mode extraction failed");
        }
        const NormalizedComplexMode normalized = normalizeAcousticMode(
            input, matrix, baseEigenvalue, raw.turningPoint, raw.mode);
        if (mode == 0)
        {
            eigen.ModeZ = raw.depth;
            eigen.PhiMode = Eigen::MatrixXcd::Zero(
                eigen.k.size(), raw.depth.size());
        }
        else if (raw.depth.size() != eigen.ModeZ.size() ||
                 !raw.depth.isApprox(eigen.ModeZ))
        {
            throw std::runtime_error("interface modes use inconsistent depth meshes");
        }
        if (normalized.mode.size() != eigen.ModeZ.size())
        {
            throw std::runtime_error("interface normalized mode depth mismatch");
        }
        eigen.PhiMode.row(mode) = normalized.mode.transpose();
        eigen.VG[mode] = normalized.groupVelocity;
        for (Eigen::Index depth = 0; depth < eigen.PsiS.cols(); ++depth)
        {
            const SampledMode sampled = sampleMode(
                raw.depth, normalized.mode,
                input.sourceDepths[static_cast<std::size_t>(depth)]);
            eigen.PsiS(mode, depth) = sampled.value;
            eigen.dPsidzS(mode, depth) = sampled.derivative;
        }
        for (Eigen::Index depth = 0; depth < eigen.PsiR.cols(); ++depth)
        {
            const SampledMode sampled = sampleMode(
                raw.depth, normalized.mode,
                input.receiverDepths[static_cast<std::size_t>(depth)]);
            eigen.PsiR(mode, depth) = sampled.value;
            eigen.dPsidzR(mode, depth) = sampled.derivative;
        }
    }
    return eigen;
}

char boundaryCode(AcousticBoundaryType type)
{
    switch (type)
    {
    case AcousticBoundaryType::Vacuum: return 'V';
    case AcousticBoundaryType::Rigid: return 'R';
    case AcousticBoundaryType::HalfSpace: return 'A';
    case AcousticBoundaryType::ReflectionCoefficient: return 'F';
    case AcousticBoundaryType::InternalReflection: return 'P';
    }
    return 'V';
}

std::vector<double> mergedModeDepths(const AcousticCase &input)
{
    std::vector<double> result = input.sourceDepths;
    result.insert(result.end(), input.receiverDepths.begin(),
                  input.receiverDepths.end());
    if (result.empty())
    {
        for (const AcousticLayer &layer : input.layers)
        {
            if (!layer.samples.empty() && layer.samples.front().cs == 0.0)
            {
                result.push_back(layer.topDepth);
                result.push_back(layer.bottomDepth);
            }
        }
    }
    std::sort(result.begin(), result.end());
    result.erase(std::unique(result.begin(), result.end(),
                             [](double left, double right) {
                                 return std::abs(left - right) <= 1.0e-8;
                             }),
                 result.end());
    return result;
}

ModeBoundaryData makeModeBoundary(const AcousticBoundary &source,
                                  double depth,
                                  const AcousticCase &input)
{
    ModeBoundaryData result;
    result.type = boundaryCode(source.type);
    result.cp = source.cp > 0.0
                    ? complexSoundSpeed(
                          depth, source.cp, source.alphaP,
                          attenuationContext(input, source.attenuationPower,
                                             source.transitionFrequency))
                    : std::complex<double>{};
    result.cs = source.cs > 0.0
                    ? complexSoundSpeed(
                          depth, source.cs, source.alphaS,
                          attenuationContext(input, source.attenuationPower,
                                             source.transitionFrequency))
                    : std::complex<double>{};
    result.rho = source.rho;
    result.depth = depth;
    return result;
}

ModeProfileData makeModeProfile(const AcousticCase &input,
                                const EigenParams &eigen)
{
    if (eigen.M < 1 || eigen.k.size() != eigen.M ||
        eigen.ModeZ.size() < 2 || eigen.PhiMode.rows() != eigen.M ||
        eigen.PhiMode.cols() != eigen.ModeZ.size())
    {
        throw std::invalid_argument("in-memory eigen result is incomplete");
    }
    ModeProfileData result;
    result.depths = mergedModeDepths(input);
    if (result.depths.size() < 2)
    {
        throw std::invalid_argument("in-memory mode sample depths are incomplete");
    }
    result.wavenumbers.resize(static_cast<std::size_t>(eigen.k.size()));
    for (Eigen::Index mode = 0; mode < eigen.k.size(); ++mode)
    {
        const std::complex<float> stored(
            static_cast<float>(eigen.k[mode].real()),
            static_cast<float>(eigen.k[mode].imag()));
        result.wavenumbers[static_cast<std::size_t>(mode)] = stored;
    }
    result.modes.resize(static_cast<std::size_t>(eigen.M));
    for (Eigen::Index mode = 0; mode < eigen.PhiMode.rows(); ++mode)
    {
        auto &target = result.modes[static_cast<std::size_t>(mode)];
        target.resize(result.depths.size());
        const Eigen::VectorXcd fullMode = eigen.PhiMode.row(mode).transpose();
        for (std::size_t depth = 0; depth < result.depths.size(); ++depth)
        {
            const std::complex<double> sampled =
                sampleMode(eigen.ModeZ, fullMode, result.depths[depth]).value;
            const std::complex<float> stored(
                static_cast<float>(sampled.real()),
                static_cast<float>(sampled.imag()));
            target[depth] = stored;
        }
    }
    std::size_t firstAcoustic = input.layers.size();
    std::size_t lastAcoustic = input.layers.size();
    for (std::size_t medium = 0; medium < input.layers.size(); ++medium)
    {
        if (!input.layers[medium].samples.empty() &&
            input.layers[medium].samples.front().cs == 0.0)
        {
            if (firstAcoustic == input.layers.size()) firstAcoustic = medium;
            lastAcoustic = medium;
        }
    }
    if (firstAcoustic == input.layers.size())
    {
        throw std::invalid_argument("in-memory modes require an acoustic layer");
    }
    for (std::size_t medium = firstAcoustic; medium <= lastAcoustic; ++medium)
    {
        const AcousticLayer &layer = input.layers[medium];
        result.meshCounts.push_back(layer.baseMesh);
        result.materials.push_back("ACOUSTIC");
        result.mediumDepths.push_back(layer.topDepth);
        result.mediumDensities.push_back(
            layer.samples.empty() ? 0.0 : layer.samples.front().rho);
    }
    result.top = makeModeBoundary(input.top,
                                  input.layers[firstAcoustic].topDepth,
                                  input);
    result.bottom = makeModeBoundary(input.bottom,
                                     input.layers[lastAcoustic].bottomDepth,
                                     input);
    return result;
}

ModeFileData makeModeData(const std::vector<AcousticCase> &inputs,
                          const std::vector<EigenParams> &eigen)
{
    if (inputs.empty() || inputs.size() != eigen.size())
    {
        throw std::invalid_argument("in-memory mode profile count mismatch");
    }
    ModeFileData result;
    result.title = inputs.front().title;
    result.frequency = inputs.front().frequency;
    static_cast<ModeProfileData &>(result) = makeModeProfile(inputs.front(), eigen.front());
    for (std::size_t profile = 1; profile < inputs.size(); ++profile)
    {
        result.additionalProfiles.push_back(makeModeProfile(inputs[profile], eigen[profile]));
    }
    return result;
}

std::string normalizedFieldParameterSignature(const OOKC_parameters &source)
{
    OOKC_parameters normalized = source;
    normalized.envPath.clear();
    normalized.flpPath.clear();
    normalized.modPath.clear();
    normalized.shdPath.clear();
    normalized.runMode = Run_Mode::MODE_B_Both;
    std::ostringstream fieldProfiles;
    fieldProfiles << "|field-profile-count=" << source.NProf << std::hexfloat;
    for (double rangeKm : source.RProf)
    {
        fieldProfiles << ':' << rangeKm;
    }
    return parameters_to_json_string(normalized) + fieldProfiles.str();
}

std::string normalizedEigenParameterSignature(const OOKC_parameters &source)
{
    OOKC_parameters normalized = source;
    normalized.envPath.clear();
    normalized.flpPath.clear();
    normalized.modPath.clear();
    normalized.shdPath.clear();
    normalized.runMode = Run_Mode::MODE_B_Both;
    normalized.Pos.Rr.resize(0);
    normalized.Pos.NRr = 0;
    normalized.Pos.is_Linspace_Rr = false;
    normalized.Pos.Ro.resize(0);
    normalized.Pos.NRo = 0;
    normalized.Pos.is_Linspace_Ro = false;
    normalized.Pos.NRz_per_range = 0;
    normalized.Pos.Delta_r = 0.0;
    normalized.Pos.GridType = Grid_Mode::MODE_R_Rectangular;
    normalized.MLimit = 1;
    normalized.SourceType = Source_Mode::MODE_R_Point;
    normalized.SBP = SrcBmPat{};
    normalized.is_Velocity = false;
    normalized.coherenceType = CoherenceType::Coherent;
    normalized.modeType = ModeType::Adiabatic;
    normalized.RProf = Eigen::VectorXd::Zero(
        static_cast<Eigen::Index>(normalized.sspInput.size()));
    for (ssp::Range_Independent_Area &profile : normalized.sspInput)
    {
        profile.Range = 0.0;
    }
    return parameters_to_json_string(normalized);
}

void clearFieldResults(OOKC_output &output)
{
    output.pressure.clear();
    output.horizontalVelocity.clear();
    output.verticalVelocity.clear();
    output.sourceCount = 0;
    output.receiverDepthCount = 0;
    output.rangeCount = 0;
}
}

class KernelInterface::InterfaceImpl
{
public:
    OOKC_parameters params;
    OOKC_output output;
    ThreadPool::ExecutorLease ownedProfileExecutor;
    ThreadPool::BorrowedExecutor borrowedProfileExecutor;
    bool hasBorrowedProfileExecutor = false;
    int fieldThreads = static_cast<int>(std::max(1u, std::thread::hardware_concurrency()));
    bool alive = true;
    std::string eigenSignature;
    std::string fieldSignature;
    std::shared_ptr<detail::ResultViewState> viewState =
        std::make_shared<detail::ResultViewState>();

    ThreadPool::ExecutorLease acquireProfileExecutor() const
    {
        if (ownedProfileExecutor)
        {
            return ownedProfileExecutor;
        }
        if (!hasBorrowedProfileExecutor)
        {
            return {};
        }
        ThreadPool::ExecutorLease lease = borrowedProfileExecutor.lock();
        if (!lease)
        {
            throw std::runtime_error("external ThreadPool expired");
        }
        return lease;
    }
};

KernelInterface::KernelInterface()
    : impl_(std::make_unique<InterfaceImpl>())
{
}

KernelInterface::KernelInterface(ThreadPool &threadPool)
    : KernelInterface()
{
    impl_->borrowedProfileExecutor = threadPool.borrowExecutor();
    impl_->hasBorrowedProfileExecutor = true;
}

KernelInterface::KernelInterface(std::shared_ptr<ThreadPool> threadPool)
    : KernelInterface()
{
    setProfileExecutor(std::move(threadPool));
}

KernelInterface::~KernelInterface() = default;
KernelInterface::KernelInterface(KernelInterface &&) noexcept = default;
KernelInterface &KernelInterface::operator=(KernelInterface &&) noexcept = default;

void KernelInterface::ensureAlive() const
{
    if (!impl_ || !impl_->alive)
    {
        throw std::runtime_error("OpenOcean-Krakenc interface has been released");
    }
}

void KernelInterface::ensureResultsFresh() const
{
    ensureAlive();
    const bool hasEigen = !impl_->output.eigen.empty();
    const bool hasField = !impl_->output.pressure.empty() ||
                          !impl_->output.horizontalVelocity.empty() ||
                          !impl_->output.verticalVelocity.empty();
    if (!hasEigen && !hasField)
    {
        return;
    }
    if ((hasEigen && impl_->eigenSignature !=
                         normalizedEigenParameterSignature(impl_->params)) ||
        (hasField && impl_->fieldSignature !=
                         normalizedFieldParameterSignature(impl_->params)))
    {
        throw std::logic_error(
            "OpenOcean-Krakenc parameters changed after the current result was computed");
    }
}

void KernelInterface::invalidateResultViews() noexcept
{
    if (impl_ && impl_->viewState)
    {
        impl_->viewState->generation.fetch_add(1, std::memory_order_acq_rel);
    }
}

void KernelInterface::setFieldThreads(int numThreads)
{
    ensureAlive();
    if (numThreads < 1)
    {
        throw std::invalid_argument("OpenOcean-Krakenc thread count must be positive");
    }
    impl_->fieldThreads = numThreads;
}

int KernelInterface::getFieldThreads() const
{
    ensureAlive();
    return impl_->fieldThreads;
}

void KernelInterface::setNumThreads(int numThreads) { setFieldThreads(numThreads); }
int KernelInterface::getNumThreads() const { return getFieldThreads(); }

int KernelInterface::getHardwareThreads() const
{
    ensureAlive();
    return static_cast<int>(std::max(1u, std::thread::hardware_concurrency()));
}

void KernelInterface::setThreadPool(ThreadPool &threadPool)
{
    ensureAlive();
    impl_->ownedProfileExecutor = {};
    impl_->borrowedProfileExecutor = threadPool.borrowExecutor();
    impl_->hasBorrowedProfileExecutor = true;
}

void KernelInterface::setProfileExecutor(std::shared_ptr<ThreadPool> threadPool)
{
    ensureAlive();
    if (!threadPool)
    {
        throw std::invalid_argument("profile ThreadPool must not be null");
    }
    impl_->ownedProfileExecutor = threadPool->acquireExecutor();
    impl_->borrowedProfileExecutor = {};
    impl_->hasBorrowedProfileExecutor = false;
}

void KernelInterface::setThreadPool(std::shared_ptr<ThreadPool> threadPool)
{
    setProfileExecutor(std::move(threadPool));
}

void KernelInterface::detachThreadPool()
{
    ensureAlive();
    impl_->ownedProfileExecutor = {};
    impl_->borrowedProfileExecutor = {};
    impl_->hasBorrowedProfileExecutor = false;
}

void KernelInterface::runEigen()
{
    ensureAlive();
    ThreadPool::ExecutorLease executor = impl_->acquireProfileExecutor();
    OOKC_parameters candidate = impl_->params;
    std::vector<AcousticCase> inputs;
    if (candidate.sspInput.empty())
    {
        if (candidate.envPath.empty())
        {
            throw std::invalid_argument(
                "OpenOcean-Krakenc eigen ENV path is not configured");
        }
        inputs = readAcousticEnvironments(candidate.envPath);
        updatePublicParameters(inputs, candidate);
    }
    else
    {
        inputs = toAcousticCases(candidate);
    }
    std::vector<EigenParams> completed;
    completed.reserve(inputs.size());
    if (executor && !executor.ownsCurrentThread())
    {
        std::vector<std::future<EigenParams>> futures;
        futures.reserve(inputs.size());
        for (const AcousticCase &input : inputs)
        {
            futures.push_back(executor.enqueue([input]() {
                return solveInterfaceProfile(input);
            }));
        }
        for (auto &future : futures)
        {
            completed.push_back(future.get());
        }
    }
    else
    {
        for (const AcousticCase &input : inputs)
        {
            completed.push_back(solveInterfaceProfile(input));
        }
    }
    invalidateResultViews();
    impl_->params = std::move(candidate);
    impl_->output.eigen = std::move(completed);
    impl_->output.pressure.clear();
    impl_->output.horizontalVelocity.clear();
    impl_->output.verticalVelocity.clear();
    impl_->output.sourceCount = 0;
    impl_->output.receiverDepthCount = 0;
    impl_->output.rangeCount = 0;
    impl_->eigenSignature = normalizedEigenParameterSignature(impl_->params);
    impl_->fieldSignature.clear();
}

void KernelInterface::runField()
{
    ensureAlive();
    const OOKC_parameters previousParams = impl_->params;
    const OOKC_output previousOutput = impl_->output;
    const std::string previousEigenSignature = impl_->eigenSignature;
    const std::string previousFieldSignature = impl_->fieldSignature;
    try
    {
    const bool hasCurrentEigen = !impl_->output.eigen.empty() &&
        impl_->eigenSignature == normalizedEigenParameterSignature(impl_->params);
    bool useEigenResults = hasCurrentEigen;
    std::filesystem::path explicitModePath = impl_->params.modPath;
    const bool hasExplicitMode = !explicitModePath.empty() &&
                                 std::filesystem::is_regular_file(explicitModePath);
    if (!useEigenResults && !hasExplicitMode)
    {
        runEigen();
        useEigenResults = true;
    }
    OOKC_parameters candidate = impl_->params;
    if (candidate.Pos.Rr.size() == 0)
    {
        std::filesystem::path flpPath = candidate.flpPath;
        if (flpPath.empty() && !candidate.envPath.empty())
        {
            flpPath = candidate.envPath;
            flpPath.replace_extension(".flp");
        }
        if (flpPath.empty() || !read_flp_file(flpPath.string(), candidate))
        {
            throw std::invalid_argument(
                "OpenOcean-Krakenc field input requires receiver ranges or a valid FLP path");
        }
    }
    std::vector<AcousticCase> inputs;
    ModeFileData modes;
    if (useEigenResults)
    {
        inputs = toAcousticCases(candidate);
        modes = makeModeData(inputs, impl_->output.eigen);
    }
    else
    {
        modes = readModeFile(explicitModePath);
        if (!candidate.sspInput.empty())
        {
            inputs = toAcousticCases(candidate);
        }
    }
    if (!(candidate.freqinfo.freq > 0.0))
    {
        candidate.freqinfo.freq = modes.frequency;
        candidate.freqinfo.Nfreq = 1;
        candidate.freqinfo.freqvec = Eigen::VectorXd::Constant(1, modes.frequency);
    }
    FieldParameters parameters = toFieldParameters(candidate);
    parameters.threadCount = impl_->fieldThreads;
    if (!inputs.empty())
    {
        for (double depth : parameters.sourceDepths)
        {
            parameters.sourceSoundSpeeds.push_back(soundSpeedAt(inputs.front(), depth));
        }
    }
    PressureField field;
    field = evaluateField(modes, parameters);
    std::vector<std::complex<float>> pressure(field.values.size());
    std::vector<std::complex<float>> horizontal(field.horizontalValues.size());
    std::vector<std::complex<float>> vertical(field.verticalValues.size());
    std::transform(field.values.begin(), field.values.end(),
                   pressure.begin(), [](std::complex<double> value) {
                       return std::complex<float>(static_cast<float>(value.real()),
                                                  static_cast<float>(value.imag()));
                   });
    std::transform(field.horizontalValues.begin(), field.horizontalValues.end(),
                   horizontal.begin(), [](std::complex<double> value) {
                       return std::complex<float>(static_cast<float>(value.real()),
                                                  static_cast<float>(value.imag()));
                   });
    std::transform(field.verticalValues.begin(), field.verticalValues.end(),
                   vertical.begin(), [](std::complex<double> value) {
                       return std::complex<float>(static_cast<float>(value.real()),
                                                  static_cast<float>(value.imag()));
                   });
    invalidateResultViews();
    impl_->params = std::move(candidate);
    impl_->output.pressure = std::move(pressure);
    if (impl_->params.is_Velocity)
    {
        impl_->output.horizontalVelocity = std::move(horizontal);
        impl_->output.verticalVelocity = std::move(vertical);
    }
    else
    {
        impl_->output.horizontalVelocity.clear();
        impl_->output.verticalVelocity.clear();
    }
    impl_->output.sourceCount = field.sourceDepths.size();
    impl_->output.receiverDepthCount = field.receiverDepths.size();
    impl_->output.rangeCount = field.rangesMetres.size();
    const std::string completedEigenSignature =
        normalizedEigenParameterSignature(impl_->params);
    const std::string completedFieldSignature =
        normalizedFieldParameterSignature(impl_->params);
    if (useEigenResults)
    {
        impl_->eigenSignature = completedEigenSignature;
    }
    else
    {
        impl_->output.eigen.clear();
        impl_->eigenSignature.clear();
    }
    impl_->fieldSignature = completedFieldSignature;
    if (!impl_->params.shdPath.empty())
    {
        exportShadeResult(impl_->params, impl_->output,
                          impl_->params.shdPath, ShadeDataType::Pressure);
    }
    }
    catch (...)
    {
        impl_->params = previousParams;
        impl_->output = previousOutput;
        impl_->eigenSignature = previousEigenSignature;
        impl_->fieldSignature = previousFieldSignature;
        throw;
    }
}

void KernelInterface::run()
{
    ensureAlive();
    const OOKC_parameters previousParams = impl_->params;
    const OOKC_output previousOutput = impl_->output;
    const std::string previousEigenSignature = impl_->eigenSignature;
    const std::string previousFieldSignature = impl_->fieldSignature;
    try
    {
        switch (impl_->params.runMode)
        {
        case Run_Mode::MODE_M_Modes:
            runEigen();
            break;
        case Run_Mode::MODE_F_Field:
            runField();
            break;
        case Run_Mode::MODE_B_Both:
            runEigen();
            runField();
            break;
        }
    }
    catch (...)
    {
        impl_->params = previousParams;
        impl_->output = previousOutput;
        impl_->eigenSignature = previousEigenSignature;
        impl_->fieldSignature = previousFieldSignature;
        throw;
    }
}

void KernelInterface::clearResults()
{
    ensureAlive();
    invalidateResultViews();
    impl_->output.clear();
    impl_->eigenSignature.clear();
    impl_->fieldSignature.clear();
}

void KernelInterface::close()
{
    if (!impl_ || !impl_->alive)
    {
        return;
    }

    invalidateResultViews();
    impl_->output.clear();
    impl_->eigenSignature.clear();
    impl_->fieldSignature.clear();
    impl_->ownedProfileExecutor = {};
    impl_->borrowedProfileExecutor = {};
    impl_->hasBorrowedProfileExecutor = false;
    impl_->alive = false;
}

bool KernelInterface::isClosed() const noexcept
{
    return !impl_ || !impl_->alive;
}

void KernelInterface::free() { close(); }

namespace
{
void assignPosition(Eigen::VectorXd &target, int &count, bool &linspace,
                    const Eigen::VectorXd &values, const char *name)
{
    if (values.size() < 1 || !values.allFinite())
    {
        throw std::invalid_argument(std::string(name) + " must be nonempty and finite");
    }
    target = values;
    count = static_cast<int>(values.size());
    linspace = false;
}

Eigen::VectorXd makeLinspace(double start, double end, int count, const char *name)
{
    if (count < 1 || !std::isfinite(start) || !std::isfinite(end))
    {
        throw std::invalid_argument(
            std::string(name) + " endpoints must be finite and count positive");
    }
    if (count > 1 && end < start)
    {
        throw std::invalid_argument(std::string(name) + " end must not precede start");
    }
    return Eigen::VectorXd::LinSpaced(count, start, end);
}

std::complex<float> *sourceView(std::vector<std::complex<float>> &values,
                                const OOKC_output &shape, int sourceIndex)
{
    const std::size_t block = shape.receiverDepthCount * shape.rangeCount;
    if (block == 0 || values.size() != shape.sourceCount * block)
    {
        throw std::logic_error("OpenOcean-Krakenc field result shape is unavailable");
    }
    if (sourceIndex < 0 || static_cast<std::size_t>(sourceIndex) >= shape.sourceCount)
    {
        throw std::out_of_range("OpenOcean-Krakenc source index is out of range");
    }
    return values.data() + static_cast<std::size_t>(sourceIndex) * block;
}

std::complex<float> *allSourcesView(
    std::vector<std::complex<float>> &values,
    const OOKC_output &shape)
{
    const std::size_t expected = shape.sourceCount *
                                 shape.receiverDepthCount * shape.rangeCount;
    if (expected == 0 || values.size() != expected)
    {
        throw std::logic_error("OpenOcean-Krakenc field result is unavailable");
    }
    return values.data();
}

const std::complex<float> *sourceView(
    const std::vector<std::complex<float>> &values,
    const OOKC_output &shape, int sourceIndex)
{
    const std::size_t block = shape.receiverDepthCount * shape.rangeCount;
    if (block == 0 || values.size() != shape.sourceCount * block)
    {
        throw std::logic_error("OpenOcean-Krakenc field result shape is unavailable");
    }
    if (sourceIndex < 0 ||
        static_cast<std::size_t>(sourceIndex) >= shape.sourceCount)
    {
        throw std::out_of_range("OpenOcean-Krakenc source index is out of range");
    }
    return values.data() + static_cast<std::size_t>(sourceIndex) * block;
}

const std::complex<float> *allSourcesView(
    const std::vector<std::complex<float>> &values,
    const OOKC_output &shape)
{
    const std::size_t expected = shape.sourceCount *
                                 shape.receiverDepthCount * shape.rangeCount;
    if (expected == 0 || values.size() != expected)
    {
        throw std::logic_error("OpenOcean-Krakenc field result is unavailable");
    }
    return values.data();
}
}

void KernelInterface::set_Title(const std::string &title)
{
    ensureAlive();
    impl_->params.Title = title;
    invalidateResultViews();
    impl_->output.clear();
}

void KernelInterface::set_Freq(double freq)
{
    ensureAlive();
    if (!std::isfinite(freq) || !(freq > 0.0))
    {
        throw std::invalid_argument("OpenOcean-Krakenc frequency must be positive");
    }
    impl_->params.freqinfo.freq = freq;
    impl_->params.freqinfo.Nfreq = 1;
    impl_->params.freqinfo.freqvec = Eigen::VectorXd::Constant(1, freq);
    invalidateResultViews();
    impl_->output.clear();
}

void KernelInterface::set_freqvec(const Eigen::VectorXd &freqvec)
{
    ensureAlive();
    if (freqvec.size() < 1 || !freqvec.allFinite() ||
        (freqvec.array() <= 0.0).any())
    {
        throw std::invalid_argument("OpenOcean-Krakenc frequencies must be positive");
    }
    impl_->params.freqinfo.freqvec = freqvec;
    impl_->params.freqinfo.Nfreq = static_cast<int>(freqvec.size());
    impl_->params.freqinfo.freq = freqvec[0];
    invalidateResultViews();
    impl_->output.clear();
}

void KernelInterface::set_SSP(const std::vector<ssp::Range_Independent_Area> &sspInput)
{
    ensureAlive();
    if (sspInput.empty())
    {
        throw std::invalid_argument("OpenOcean-Krakenc SSP input must not be empty");
    }
    OOKC_parameters candidate = impl_->params;
    candidate.sspInput = sspInput;
    candidate.NProf = static_cast<int>(sspInput.size());
    candidate.RProf.resize(candidate.NProf);
    for (Eigen::Index i = 0; i < candidate.RProf.size(); ++i)
    {
        candidate.RProf[i] = sspInput[static_cast<std::size_t>(i)].Range;
    }
    candidate.freqinfo.freq = 1.0;
    candidate.freqinfo.Nfreq = 1;
    candidate.freqinfo.freqvec = Eigen::VectorXd::Constant(1, 1.0);
    candidate.AttenUnit.absModel = OceanAbsorptionModel::None;
    candidate.cLow = 1.0;
    candidate.cHigh = 2.0;
    candidate.hasModePos = false;
    candidate.Pos.Sz = Eigen::VectorXd::Zero(1);
    candidate.Pos.NSz = 1;
    candidate.Pos.Rz = Eigen::VectorXd::Zero(1);
    candidate.Pos.NRz = 1;
    validatePublicParameters(candidate, Run_Mode::MODE_M_Modes);
    impl_->params.sspInput = std::move(candidate.sspInput);
    impl_->params.NProf = candidate.NProf;
    impl_->params.RProf = std::move(candidate.RProf);
    invalidateResultViews();
    impl_->output.clear();
}

void KernelInterface::set_AttenUnit(Atten_Mode mode) { ensureAlive(); impl_->params.AttenUnit = mode; invalidateResultViews(); impl_->output.clear(); }
void KernelInterface::set_Sz(const Eigen::VectorXd &v) { ensureAlive(); assignPosition(impl_->params.Pos.Sz, impl_->params.Pos.NSz, impl_->params.Pos.is_Linspace_Sz, v, "Sz"); invalidateResultViews(); impl_->output.clear(); }
void KernelInterface::set_Rr(const Eigen::VectorXd &v) { ensureAlive(); assignPosition(impl_->params.Pos.Rr, impl_->params.Pos.NRr, impl_->params.Pos.is_Linspace_Rr, v, "Rr"); invalidateResultViews(); clearFieldResults(impl_->output); impl_->fieldSignature.clear(); }
void KernelInterface::set_Rz(const Eigen::VectorXd &v) { ensureAlive(); assignPosition(impl_->params.Pos.Rz, impl_->params.Pos.NRz, impl_->params.Pos.is_Linspace_Rz, v, "Rz"); invalidateResultViews(); impl_->output.clear(); }
void KernelInterface::set_Ro(const Eigen::VectorXd &v) { ensureAlive(); assignPosition(impl_->params.Pos.Ro, impl_->params.Pos.NRo, impl_->params.Pos.is_Linspace_Ro, v, "Ro"); invalidateResultViews(); clearFieldResults(impl_->output); impl_->fieldSignature.clear(); }
void KernelInterface::set_Sz(double a, double b, int n) { set_Sz(makeLinspace(a, b, n, "Sz")); impl_->params.Pos.is_Linspace_Sz = true; }
void KernelInterface::set_Rr(double a, double b, int n) { set_Rr(makeLinspace(a, b, n, "Rr")); impl_->params.Pos.is_Linspace_Rr = true; }
void KernelInterface::set_Rz(double a, double b, int n) { set_Rz(makeLinspace(a, b, n, "Rz")); impl_->params.Pos.is_Linspace_Rz = true; }
void KernelInterface::set_Ro(double a, double b, int n) { set_Ro(makeLinspace(a, b, n, "Ro")); impl_->params.Pos.is_Linspace_Ro = true; }

void KernelInterface::set_cPhase(double cLow, double cHigh)
{
    ensureAlive();
    if (!std::isfinite(cLow) || !std::isfinite(cHigh) ||
        cLow < 0.0 || !(cHigh > cLow))
    {
        throw std::invalid_argument("OpenOcean-Krakenc phase speed limits are invalid");
    }
    impl_->params.cLow = cLow;
    impl_->params.cHigh = cHigh;
    invalidateResultViews();
    impl_->output.clear();
}

void KernelInterface::set_GridType(Grid_Mode type) { ensureAlive(); impl_->params.Pos.GridType = type; invalidateResultViews(); clearFieldResults(impl_->output); impl_->fieldSignature.clear(); }
void KernelInterface::set_Rmax(double rMax) { ensureAlive(); if (!std::isfinite(rMax) || rMax < 0.0) throw std::invalid_argument("Rmax must be finite and nonnegative"); impl_->params.Rmax = rMax; invalidateResultViews(); impl_->output.clear(); }
void KernelInterface::set_SourceType(Source_Mode type) { ensureAlive(); impl_->params.SourceType = type; invalidateResultViews(); clearFieldResults(impl_->output); impl_->fieldSignature.clear(); }
void KernelInterface::set_MLimit(int limit)
{
    ensureAlive();
    if (limit <= 0)
    {
        throw std::invalid_argument("mode limit must be positive");
    }
    impl_->params.MLimit = limit;
    invalidateResultViews();
    impl_->output.clear();
    impl_->eigenSignature.clear();
    impl_->fieldSignature.clear();
}
void KernelInterface::set_CoherenceType(CoherenceType type)
{
    ensureAlive();
    impl_->params.coherenceType = type;
    invalidateResultViews();
    clearFieldResults(impl_->output);
    impl_->fieldSignature.clear();
}
void KernelInterface::set_ModeType(ModeType type)
{
    ensureAlive();
    impl_->params.modeType = type;
    invalidateResultViews();
    clearFieldResults(impl_->output);
    impl_->fieldSignature.clear();
}
void KernelInterface::set_ModeSampling(
    const Eigen::VectorXd &sourceDepths,
    const Eigen::VectorXd &receiverDepths)
{
    ensureAlive();
    if (sourceDepths.size() < 1 || receiverDepths.size() < 1 ||
        !sourceDepths.allFinite() || !receiverDepths.allFinite())
    {
        throw std::invalid_argument(
            "mode sampling depths must be nonempty and finite");
    }
    Position mode;
    mode.NSz = static_cast<int>(sourceDepths.size());
    mode.NRz = static_cast<int>(receiverDepths.size());
    mode.NRz_per_range = mode.NRz;
    mode.Sz = sourceDepths;
    mode.Rz = receiverDepths;
    mode.GridType = Grid_Mode::MODE_R_Rectangular;
    impl_->params.ModePos = std::move(mode);
    impl_->params.hasModePos = true;
    invalidateResultViews();
    impl_->output.clear();
    impl_->eigenSignature.clear();
    impl_->fieldSignature.clear();
}
void KernelInterface::set_RunMode(Run_Mode mode) { ensureAlive(); impl_->params.runMode = mode; }
void KernelInterface::set_Velocity_enable(bool enabled) { ensureAlive(); impl_->params.is_Velocity = enabled; invalidateResultViews(); clearFieldResults(impl_->output); impl_->fieldSignature.clear(); }

void KernelInterface::set_ReflCoef_Top(const std::vector<ReflectionCoef> &values)
{
    ensureAlive();
    for (std::size_t i = 0; i < values.size(); ++i)
    {
        if (!std::isfinite(values[i].theta) || !std::isfinite(values[i].R) ||
            !std::isfinite(values[i].phi) || values[i].R < 0.0 ||
            (i > 0 && values[i].theta <= values[i - 1].theta))
        {
            throw std::invalid_argument(
                "top reflection coefficients must be finite and angle ordered");
        }
    }
    impl_->params.ReflectionCoef.RTop.resize(1, static_cast<Eigen::Index>(values.size()));
    for (Eigen::Index i = 0; i < impl_->params.ReflectionCoef.RTop.size(); ++i) impl_->params.ReflectionCoef.RTop[i] = values[static_cast<std::size_t>(i)];
    impl_->params.ReflectionCoef.isDeg = false;
    invalidateResultViews();
    impl_->output.clear();
}

void KernelInterface::set_ReflCoef_Bottom(const std::vector<ReflectionCoef> &values)
{
    ensureAlive();
    for (std::size_t i = 0; i < values.size(); ++i)
    {
        if (!std::isfinite(values[i].theta) || !std::isfinite(values[i].R) ||
            !std::isfinite(values[i].phi) || values[i].R < 0.0 ||
            (i > 0 && values[i].theta <= values[i - 1].theta))
        {
            throw std::invalid_argument(
                "bottom reflection coefficients must be finite and angle ordered");
        }
    }
    impl_->params.ReflectionCoef.RBot.resize(1, static_cast<Eigen::Index>(values.size()));
    for (Eigen::Index i = 0; i < impl_->params.ReflectionCoef.RBot.size(); ++i) impl_->params.ReflectionCoef.RBot[i] = values[static_cast<std::size_t>(i)];
    impl_->params.ReflectionCoef.isDeg = false;
    invalidateResultViews();
    impl_->output.clear();
}

void KernelInterface::set_SBP(const Eigen::VectorXd &pattern, const Eigen::VectorXd &anglesDegrees)
{
    ensureAlive();
    if (pattern.size() < 2 || pattern.size() != anglesDegrees.size() ||
        !pattern.allFinite() || !anglesDegrees.allFinite() ||
        (anglesDegrees.tail(anglesDegrees.size() - 1) -
         anglesDegrees.head(anglesDegrees.size() - 1)).minCoeff() <= 0.0)
    {
        throw std::invalid_argument(
            "OpenOcean-Krakenc source beam pattern must be finite and strictly ordered");
    }
    impl_->params.SBP.pat = pattern;
    impl_->params.SBP.theta = anglesDegrees;
    impl_->params.SBP.NSBPPts = static_cast<int>(pattern.size());
    impl_->params.SBP.isSet = true;
    invalidateResultViews();
    clearFieldResults(impl_->output);
    impl_->fieldSignature.clear();
}

std::complex<float> *KernelInterface::get_u(int i) { ensureResultsFresh(); return sourceView(impl_->output.pressure, impl_->output, i); }
std::complex<float> *KernelInterface::get_v(int i) { ensureResultsFresh(); return sourceView(impl_->output.verticalVelocity, impl_->output, i); }
std::complex<float> *KernelInterface::get_h(int i) { ensureResultsFresh(); return sourceView(impl_->output.horizontalVelocity, impl_->output, i); }
std::complex<float> *KernelInterface::get_u_AllSources() { ensureResultsFresh(); return allSourcesView(impl_->output.pressure, impl_->output); }
std::complex<float> *KernelInterface::get_v_AllSources() { ensureResultsFresh(); return allSourcesView(impl_->output.verticalVelocity, impl_->output); }
std::complex<float> *KernelInterface::get_h_AllSources() { ensureResultsFresh(); return allSourcesView(impl_->output.horizontalVelocity, impl_->output); }

ArrayView<const std::complex<float>> KernelInterface::get_u_view(int i) const
{
    ensureResultsFresh();
    const OOKC_output &output = impl_->output;
    const std::size_t size = output.receiverDepthCount * output.rangeCount;
    return {sourceView(output.pressure, output, i), size,
            impl_->viewState,
            impl_->viewState->generation.load(std::memory_order_acquire)};
}

ArrayView<const std::complex<float>> KernelInterface::get_v_view(int i) const
{
    ensureResultsFresh();
    const OOKC_output &output = impl_->output;
    const std::size_t size = output.receiverDepthCount * output.rangeCount;
    return {sourceView(output.verticalVelocity, output, i), size,
            impl_->viewState,
            impl_->viewState->generation.load(std::memory_order_acquire)};
}

ArrayView<const std::complex<float>> KernelInterface::get_h_view(int i) const
{
    ensureResultsFresh();
    const OOKC_output &output = impl_->output;
    const std::size_t size = output.receiverDepthCount * output.rangeCount;
    return {sourceView(output.horizontalVelocity, output, i), size,
            impl_->viewState,
            impl_->viewState->generation.load(std::memory_order_acquire)};
}

ArrayView<const std::complex<float>> KernelInterface::get_u_view_all_sources() const
{
    ensureResultsFresh();
    const OOKC_output &output = impl_->output;
    return {allSourcesView(output.pressure, output),
            output.pressure.size(), impl_->viewState,
            impl_->viewState->generation.load(std::memory_order_acquire)};
}

ArrayView<const std::complex<float>> KernelInterface::get_v_view_all_sources() const
{
    ensureResultsFresh();
    const OOKC_output &output = impl_->output;
    return {allSourcesView(output.verticalVelocity, output),
            output.verticalVelocity.size(), impl_->viewState,
            impl_->viewState->generation.load(std::memory_order_acquire)};
}

ArrayView<const std::complex<float>> KernelInterface::get_h_view_all_sources() const
{
    ensureResultsFresh();
    const OOKC_output &output = impl_->output;
    return {allSourcesView(output.horizontalVelocity, output),
            output.horizontalVelocity.size(), impl_->viewState,
            impl_->viewState->generation.load(std::memory_order_acquire)};
}

void KernelInterface::export_result(const std::string &root)
{
    ensureAlive();
    export_mod(root + "_P");
    export_shd(root + "_P", 1);
    if (impl_->params.is_Velocity)
    {
        export_shd(root + "_V", 2);
        export_shd(root + "_H", 3);
    }
}
void KernelInterface::export_mod(const std::string &root)
{
    ensureResultsFresh();
    exportModeResult(impl_->params, impl_->output, resultPath(root, ".mod"));
}
void KernelInterface::export_shd(const std::string &root, int dataType)
{
    ensureResultsFresh();
    ShadeDataType type;
    switch (dataType)
    {
    case 1: type = ShadeDataType::Pressure; break;
    case 2: type = ShadeDataType::VerticalVelocity; break;
    case 3: type = ShadeDataType::HorizontalVelocity; break;
    default: throw std::invalid_argument("SHD dataType must be 1, 2, or 3");
    }
    exportShadeResult(impl_->params, impl_->output,
                      resultPath(root, ".shd"), type);
}
LoadResult KernelInterface::loadEnv(const std::string &envPath)
{
    ensureAlive();
    OOKC_parameters candidate;
    const LoadResult envResult = read_env_file_result(envPath, candidate);
    if (!envResult.ok)
    {
        return envResult;
    }
    const LoadResult flpResult = read_flp_file_result(envPath, candidate);
    if (!flpResult.ok)
    {
        return flpResult;
    }
    try
    {
        validatePublicParameters(candidate, Run_Mode::MODE_B_Both);
    }
    catch (const std::exception &error)
    {
        return LoadResult::failure(LoadErrorCode::InvalidField, envPath,
                                   error.what());
    }
    invalidateResultViews();
    impl_->params = std::move(candidate);
    impl_->output.clear();
    impl_->eigenSignature.clear();
    impl_->fieldSignature.clear();
    return LoadResult::success();
}

bool KernelInterface::from_env(const std::string &envPath)
{
    return loadEnv(envPath).ok;
}

#if defined(OPENOCEAN_KRAKENC_LEGACY_JSON)
LoadResult KernelInterface::loadJson(const std::string &jsonPath)
{
    ensureAlive();
    OOKC_parameters candidate;
    const LoadResult result = read_json_file_result(jsonPath, candidate);
    if (!result.ok)
    {
        return result;
    }
    invalidateResultViews();
    impl_->params = std::move(candidate);
    impl_->output.clear();
    impl_->eigenSignature.clear();
    impl_->fieldSignature.clear();
    return LoadResult::success();
}

bool KernelInterface::from_json(const std::string &jsonPath)
{
    return loadJson(jsonPath).ok;
}
bool KernelInterface::to_json(const std::string &jsonPath) const { ensureAlive(); return write_json_file(jsonPath, impl_->params); }
std::string KernelInterface::to_json_string() const { ensureAlive(); return parameters_to_json_string(impl_->params); }
#endif

OOKC_parameters &KernelInterface::getParams()
{
    ensureAlive();
    invalidateResultViews();
    return impl_->params;
}

const OOKC_parameters &KernelInterface::getParams() const
{
    ensureAlive();
    return impl_->params;
}

const OOKC_parameters &KernelInterface::getParams_const() const { return getParams(); }

OOKC_output &KernelInterface::getOutput()
{
    ensureResultsFresh();
    invalidateResultViews();
    return impl_->output;
}

const OOKC_output &KernelInterface::getOutput() const
{
    ensureResultsFresh();
    return impl_->output;
}

const OOKC_output &KernelInterface::getOutput_const() const { return getOutput(); }
OOKC_output KernelInterface::getOutput_Copy() const { return getOutput(); }
}
