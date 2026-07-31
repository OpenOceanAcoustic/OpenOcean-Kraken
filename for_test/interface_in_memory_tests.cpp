#include "OpenOceanKrakencKernelInterface.h"
#include "module/env_in_out.hpp"
#include "algorithm/FieldSolver.h"

#include <cassert>
#include <complex>
#include <filesystem>

using namespace OpenOceanKrakenc;

#undef assert
#define assert(condition) do { if (!(condition)) throw std::runtime_error("assertion failed: " #condition); } while (false)

namespace
{
void require(bool condition, const char *message)
{
    if (!condition) throw std::runtime_error(message);
}
}

int main()
{
    const std::filesystem::path root = OPENOCEANKRAKENC_SOURCE_DIR;
    const std::filesystem::path env = root / "for_test" / "fixtures" / "two_profile_small.env";
    const std::filesystem::path flp = root / "for_test" / "fixtures" / "two_profile_small.flp";
    const FieldParameters parsed = readFieldParameters(flp);
    require(parsed.profileRangesKm.size() == 2, "test FLP profile count mismatch");

    KernelInterface memory;
    require(memory.from_env(env.string()), "from_env failed for in-memory fixture");
    memory.set_Velocity_enable(true);
    require(read_flp_file(flp.string(), memory.getParams()), "read_flp_file failed for in-memory fixture");
    memory.getParams().envPath.clear();
    memory.getParams().flpPath.clear();
    memory.getParams().modPath.clear();
    memory.getParams().shdPath.clear();
    memory.run();

    const OOKC_output &memoryOutput = memory.getOutput_const();
    assert(memoryOutput.eigen.size() == 2);
    assert(memoryOutput.eigen.front().M > 0);
    assert(memoryOutput.eigen.front().ModeZ.size() >= 2);
    assert(memoryOutput.eigen.front().PhiMode.rows() == memoryOutput.eigen.front().M);
    assert(memoryOutput.eigen.front().PhiMode.cols() == memoryOutput.eigen.front().ModeZ.size());
    assert(memoryOutput.sourceCount == 1);
    assert(memoryOutput.receiverDepthCount == 1);
    assert(memoryOutput.rangeCount == 1);
    assert(memoryOutput.pressure.size() == 1);
    assert(memory.get_u(0) == memory.get_u_AllSources());
    assert(memory.get_v(0) == memory.get_v_AllSources());
    assert(memory.get_h(0) == memory.get_h_AllSources());

    const std::filesystem::path jsonPath =
        std::filesystem::temp_directory_path() /
        "openocean_krakenc_interface_equivalence.json";
    require(memory.to_json(jsonPath.string()), "in-memory JSON write failed");
    KernelInterface json;
    require(json.from_json(jsonPath.string()), "in-memory JSON load failed");
    json.run();
    const OOKC_output &jsonOutput = json.getOutput_const();
    assert(jsonOutput.eigen.size() == memoryOutput.eigen.size());
    for (std::size_t profile = 0; profile < memoryOutput.eigen.size(); ++profile)
    {
        assert(jsonOutput.eigen[profile].k == memoryOutput.eigen[profile].k);
        assert(jsonOutput.eigen[profile].VG == memoryOutput.eigen[profile].VG);
        assert(jsonOutput.eigen[profile].PhiMode ==
               memoryOutput.eigen[profile].PhiMode);
    }
    assert(jsonOutput.pressure == memoryOutput.pressure);
    assert(jsonOutput.verticalVelocity == memoryOutput.verticalVelocity);
    assert(jsonOutput.horizontalVelocity == memoryOutput.horizontalVelocity);

    const std::filesystem::path modeRoot =
        std::filesystem::temp_directory_path() /
        "openocean_krakenc_explicit_mode_input";
    const std::filesystem::path shadePath =
        std::filesystem::temp_directory_path() /
        "openocean_krakenc_automatic_path_output.shd";
    memory.export_mod(modeRoot.string());
    KernelInterface modeFileInput;
    modeFileInput.getParams().modPath = modeRoot.string() + ".mod";
    modeFileInput.getParams().flpPath = flp.string();
    modeFileInput.getParams().shdPath = shadePath.string();
    modeFileInput.runField();
    assert(modeFileInput.getOutput_const().eigen.empty());
    assert(!modeFileInput.getOutput_const().pressure.empty());
    assert(std::abs(modeFileInput.getParams_const().freqinfo.freq - 50.0) < 1.0e-12);
    assert(std::filesystem::exists(shadePath));
    modeFileInput.getParams().RProf[0] = 1.0;
    bool staleProfileRangeRejected = false;
    try
    {
        static_cast<void>(modeFileInput.getOutput_const());
    }
    catch (const std::logic_error &)
    {
        staleProfileRangeRejected = true;
    }
    assert(staleProfileRangeRejected);
    modeFileInput.getParams().RProf[0] = 0.0;

    KernelInterface legacy;
    legacy.getParams().envPath = env.string();
    legacy.getParams().flpPath = flp.string();
    legacy.run();
    const OOKC_output &legacyOutput = legacy.getOutput_const();
    assert(legacyOutput.eigen.size() == memoryOutput.eigen.size());
    assert(legacyOutput.pressure.size() == memoryOutput.pressure.size());
    assert((legacyOutput.eigen.front().k - memoryOutput.eigen.front().k).norm() < 1.0e-12);
    assert(std::abs(legacyOutput.pressure.front() - memoryOutput.pressure.front()) < 1.0e-6f);

    KernelInterface emptyResults;
    bool emptyAllSourcesRejected = false;
    try
    {
        static_cast<void>(emptyResults.get_u_AllSources());
    }
    catch (const std::logic_error &)
    {
        emptyAllSourcesRejected = true;
    }
    assert(emptyAllSourcesRejected);

    KernelInterface pressureOnly;
    require(pressureOnly.from_env(env.string()), "pressure-only from_env failed");
    pressureOnly.getParams().shdPath.clear();
    pressureOnly.runField();
    assert(!pressureOnly.getOutput_const().pressure.empty());
    assert(pressureOnly.getOutput_const().verticalVelocity.empty());
    assert(pressureOnly.getOutput_const().horizontalVelocity.empty());
    bool unavailableVelocityRejected = false;
    try
    {
        static_cast<void>(pressureOnly.get_v_AllSources());
    }
    catch (const std::logic_error &)
    {
        unavailableVelocityRejected = true;
    }
    assert(unavailableVelocityRejected);
    const std::size_t retainedEigenProfiles =
        pressureOnly.getOutput_const().eigen.size();
    pressureOnly.set_Rr(Eigen::VectorXd::Constant(1, 1500.0));
    assert(pressureOnly.getOutput_const().eigen.size() == retainedEigenProfiles);
    assert(pressureOnly.getOutput_const().pressure.empty());

    KernelInterface belowSsp;
    require(belowSsp.from_env(env.string()), "below-SSP from_env failed");
    belowSsp.set_Sz(Eigen::VectorXd::Constant(1, 150.0));
    belowSsp.set_Rz(Eigen::VectorXd::Constant(1, 150.0));
    belowSsp.getParams().hasModePos = false;
    belowSsp.runEigen();
    assert(belowSsp.getOutput_const().eigen.front().PsiS.cwiseAbs().maxCoeff() == 0.0);
    assert(belowSsp.getOutput_const().eigen.front().PsiR.cwiseAbs().maxCoeff() == 0.0);

    KernelInterface stale;
    require(stale.from_env(env.string()), "stale-result from_env failed");
    stale.getParams().shdPath.clear();
    stale.run();
    const Eigen::VectorXcd oldWavenumbers = stale.getOutput_const().eigen.front().k;
    stale.getParams().freqinfo.freq = 60.0;
    stale.getParams().freqinfo.freqvec = Eigen::VectorXd::Constant(1, 60.0);
    bool staleRejected = false;
    try
    {
        static_cast<void>(stale.getOutput_const());
    }
    catch (const std::logic_error &)
    {
        staleRejected = true;
    }
    assert(staleRejected);
    stale.runField();
    const Eigen::VectorXcd &newWavenumbers =
        stale.getOutput_const().eigen.front().k;
    assert(newWavenumbers.size() != oldWavenumbers.size() ||
           (newWavenumbers - oldWavenumbers).norm() > 1.0e-8);

    KernelInterface transactional;
    require(transactional.from_env(env.string()), "transactional from_env failed");
    transactional.set_GridType(Grid_Mode::MODE_I_Irregular);
    bool rejected = false;
    try
    {
        transactional.runField();
    }
    catch (const std::invalid_argument &)
    {
        rejected = true;
    }
    assert(rejected);
    assert(transactional.getOutput_const().eigen.empty());
    assert(transactional.getOutput_const().pressure.empty());

    const OOKC_output copy = memory.getOutput_Copy();
    memory.set_Title("dirty state");
    assert(memory.getOutput_const().pressure.empty());
    assert(!copy.pressure.empty());
    std::error_code ignored;
    std::filesystem::remove(jsonPath, ignored);
    std::filesystem::remove(modeRoot.string() + ".mod", ignored);
    std::filesystem::remove(shadePath, ignored);
    return 0;
}
