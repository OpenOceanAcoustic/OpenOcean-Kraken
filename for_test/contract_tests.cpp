#include "OpenOceanKrakencParams.h"
#include "OpenOceanKrakencKernelInterface.h"
#include "ThreadPool.h"

#include <algorithm>
#include <chrono>
#include <complex>
#include <future>
#include <iostream>
#include <stdexcept>
#include <string>
#include <type_traits>

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
    static_assert(std::is_same_v<EigenParams::ComplexScalar, std::complex<double>>);
    static_assert(std::is_same_v<TridMtx::ComplexScalar, std::complex<double>>);

    EigenParams eigen;
    eigen.resize(12, 5, 2, 3);
    require(eigen.firstM == 12, "modal capacity mismatch");
    require(eigen.EVMat.size() == 60, "EVMat must hold every mesh set");
    require(eigen.Extrap.size() == 60, "Extrap must hold every mesh set");
    require(eigen.k.size() == 12, "wavenumber capacity mismatch");
    require(eigen.PsiS.rows() == 12 && eigen.PsiS.cols() == 2, "source mode shape mismatch");
    require(eigen.PsiR.rows() == 12 && eigen.PsiR.cols() == 3, "receiver mode shape mismatch");

    TridMtx trid;
    trid.resize(101, 2, 5);
    require(trid.B1.size() == 101, "B1 size mismatch");
    require(trid.B4.size() == 101, "B4 size mismatch");
    require(trid.rho.size() == 101, "rho size mismatch");
    require(trid.N.size() == 2 && trid.hV.size() == 5, "mesh metadata mismatch");

    bool rejected = false;
    try
    {
        eigen.resize(0, 5, 1, 1);
    }
    catch (const std::invalid_argument &)
    {
        rejected = true;
    }
    require(rejected, "zero modal capacity must be rejected");

    ThreadPool pool(2);
    KernelInterface api(pool);
    require(api.getNumThreads() >= 1, "default thread count must be positive");
    api.setNumThreads(2);
    require(api.getNumThreads() == 2, "thread count setter failed");

    bool missingInput = false;
    try
    {
        api.runEigen();
    }
    catch (const std::invalid_argument &error)
    {
        missingInput = std::string(error.what()) ==
                       "OpenOcean-Krakenc eigen ENV path is not configured";
    }
    require(missingInput, "runEigen must reject a missing ENV path");

    api.getParams().envPath =
        std::string(OPENOCEANKRAKENC_WORKSPACE_DIR) +
        "/OpenOcean-Krakenc/for_test/fixtures/two_profile_small.env";
    api.runEigen();
    require(api.getOutput().eigen.size() == 2,
            "multi-profile interface eigen output count mismatch");
    require(std::all_of(api.getOutput().eigen.begin(),
                        api.getOutput().eigen.end(), [](const EigenParams &profile) {
                            return profile.M > 0 && profile.k.size() == profile.M &&
                                   profile.VG.size() == profile.M &&
                                   profile.VG.allFinite() && profile.VG.minCoeff() > 0.0 &&
                                   profile.PsiS.rows() == profile.M &&
                                   profile.PsiS.cols() == 1 &&
                                   profile.PsiR.rows() == profile.M &&
                                   profile.PsiR.cols() == 3 &&
                                   profile.PsiR.norm() > 0.0 &&
                                   profile.dPsidzS.rows() == profile.M &&
                                   profile.dPsidzR.cols() == 3 &&
                                   profile.dPsidzR.norm() > 0.0;
                        }),
            "multi-profile interface eigen output is incomplete");

    ThreadPool singleWorkerPool(1);
    KernelInterface reentrant(singleWorkerPool);
    reentrant.getParams().envPath =
        std::string(OPENOCEANKRAKENC_WORKSPACE_DIR) +
        "/OpenOcean-Krakenc/for_test/fixtures/two_profile_small.env";
    auto reentrantRun = singleWorkerPool.enqueue([&reentrant]() {
        reentrant.runEigen();
        return reentrant.getOutput().eigen.size();
    });
    require(reentrantRun.wait_for(std::chrono::seconds(10)) ==
                std::future_status::ready && reentrantRun.get() == 2,
            "single-worker injected pool deadlocked during reentrant runEigen");

    ThreadPool identifiedPool(1);
    auto identifiedTask = identifiedPool.enqueue_with_id("contract", []() {
        return 42;
    });
    require(identifiedTask.get() == 42,
            "identified thread-pool task returned the wrong result");
    require(identifiedPool.wait_id_for("contract", std::chrono::milliseconds(100)),
            "identified thread-pool task did not clear its ID counter");
    require(identifiedPool.wait_completion_for(std::chrono::milliseconds(100)),
            "identified thread-pool task corrupted the global completion counter");

    api.clearResults();
    api.free();
    bool released = false;
    try
    {
        static_cast<void>(api.getParams());
    }
    catch (const std::runtime_error &)
    {
        released = true;
    }
    require(released, "released interface must reject access");

    std::cout << "CONTRACT_TESTS_OK\n";
}
