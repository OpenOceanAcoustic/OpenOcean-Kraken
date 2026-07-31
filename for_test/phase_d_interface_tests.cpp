#include "OpenOceanKrakencKernelInterface.h"
#include "ThreadPool.h"

#include <chrono>
#include <complex>
#include <filesystem>
#include <fstream>
#include <future>
#include <memory>
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

std::filesystem::path fixture(const char *name)
{
    return std::filesystem::path(OPENOCEANKRAKENC_SOURCE_DIR) /
           "for_test" / "fixtures" / name;
}
}

int main()
{
    static_assert(std::is_constructible_v<KernelInterface, std::shared_ptr<ThreadPool>>);
    static_assert(std::is_same_v<decltype(std::declval<KernelInterface &>().get_u_view(0)),
                                 ArrayView<const std::complex<float>>>);
    static_assert(std::is_same_v<decltype(std::declval<KernelInterface &>().loadJson("")),
                                 LoadResult>);

    const std::filesystem::path env = fixture("two_profile_small.env");

    auto ownedPool = std::make_shared<ThreadPool>(2);
    KernelInterface owned(ownedPool);
    require(owned.loadEnv(env.string()).ok, "owning interface ENV load failed");
    ownedPool.reset();
    owned.runEigen();
    require(owned.getOutput_const().eigen.size() == 2,
            "owning pool API did not retain the executor");

    auto survivingPool = std::make_shared<ThreadPool>(1);
    {
        KernelInterface shorterLivedInterface(survivingPool);
    }
    require(survivingPool->enqueue([] { return 23; }).get() == 23,
            "destroying an interface stopped the caller-owned ThreadPool");

    auto borrowedPool = std::make_unique<ThreadPool>(2);
    KernelInterface borrowed(*borrowedPool);
    require(borrowed.loadEnv(env.string()).ok, "borrowed interface ENV load failed");
    borrowedPool.reset();
    bool expiredRejected = false;
    try
    {
        borrowed.runEigen();
    }
    catch (const std::runtime_error &error)
    {
        expiredRejected = std::string(error.what()) == "external ThreadPool expired";
    }
    require(expiredRejected, "expired borrowed pool was not rejected safely");

    auto reboundPool = std::make_shared<ThreadPool>(2);
    borrowed.setProfileExecutor(reboundPool);
    borrowed.runEigen();
    borrowed.detachThreadPool();
    reboundPool.reset();
    borrowed.runEigen();
    require(borrowed.getOutput_const().eigen.size() == 2,
            "detached interface did not fall back to serial profile solves");

    auto leaseOwner = std::make_unique<ThreadPool>(1);
    auto lease = leaseOwner->acquireExecutor();
    leaseOwner.reset();
    require(lease.enqueue([] { return 17; }).get() == 17,
            "executor lease did not survive ThreadPool object destruction");

    KernelInterface safe;
    const LoadResult envLoad = safe.loadEnv(env.string());
    require(envLoad.ok, "structured ENV load failed");
    safe.set_Velocity_enable(true);
    safe.setFieldThreads(1);
    require(safe.getFieldThreads() == 1, "field thread count was not observable");
    safe.run();
    const ArrayView<const std::complex<float>> pressure = safe.get_u_view(0);
    require(pressure.size() == 1 && pressure.at(0) == safe.getOutput_const().pressure[0],
            "safe pressure view has the wrong extent or value");

    const OOKC_output oneThread = safe.getOutput_Copy();
    safe.setFieldThreads(4);
    safe.runField();
    require(safe.getFieldThreads() == 4, "four-thread setting was not retained");
    require(safe.getOutput_const().pressure == oneThread.pressure &&
                safe.getOutput_const().verticalVelocity == oneThread.verticalVelocity &&
                safe.getOutput_const().horizontalVelocity == oneThread.horizontalVelocity,
            "field output changed with the configured thread count");

    safe.set_Rr(Eigen::VectorXd::Constant(1, 1500.0));
    bool invalidated = false;
    try
    {
        static_cast<void>(pressure.at(0));
    }
    catch (const std::logic_error &)
    {
        invalidated = true;
    }
    require(invalidated, "result view was not invalidated after result mutation");

    ArrayView<const std::complex<float>> destroyedOwnerView;
    {
        KernelInterface temporary;
        require(temporary.loadEnv(env.string()).ok,
                "temporary interface ENV load failed");
        temporary.runField();
        destroyedOwnerView = temporary.get_u_view_all_sources();
    }
    bool destroyedOwnerRejected = false;
    try
    {
        static_cast<void>(destroyedOwnerView.at(0));
    }
    catch (const std::logic_error &)
    {
        destroyedOwnerRejected = true;
    }
    require(destroyedOwnerRejected,
            "result view outlived its owning interface without rejection");

    const std::filesystem::path temp = std::filesystem::temp_directory_path();
    const LoadResult missing = safe.loadJson((temp / "ookc_missing_phase_d.json").string());
    require(!missing.ok && missing.error.code == LoadErrorCode::MissingFile,
            "missing JSON file was not classified");

    const std::filesystem::path malformed = temp / "ookc_malformed_phase_d.json";
    {
        std::ofstream stream(malformed);
        stream << "{broken";
    }
    const LoadResult malformedResult = safe.loadJson(malformed.string());
    require(!malformedResult.ok && malformedResult.error.code == LoadErrorCode::ParseError &&
                !malformedResult.error.message.empty(),
            "malformed JSON error details were not preserved");

    const std::filesystem::path unsupported = temp / "ookc_unsupported_phase_d.json";
    {
        std::ofstream stream(unsupported);
        stream << R"({"schema":"not-ookc","schema_version":1})";
    }
    const LoadResult unsupportedResult = safe.loadJson(unsupported.string());
    require(!unsupportedResult.ok &&
                unsupportedResult.error.code == LoadErrorCode::UnsupportedCapability,
            "unsupported schema was not classified");

    const std::filesystem::path incomplete = temp / "ookc_incomplete_phase_d.json";
    {
        std::ofstream stream(incomplete);
        stream << R"({"schema":"OpenOcean-Krakenc.parameters","schema_version":1})";
    }
    const LoadResult incompleteResult = safe.loadJson(incomplete.string());
    require(!incompleteResult.ok && incompleteResult.error.code == LoadErrorCode::SchemaError,
            "missing required JSON field was not classified");

    safe.close();
    safe.close();
    require(safe.isClosed(), "close must be idempotent and observable");
    bool closedRejected = false;
    try
    {
        safe.run();
    }
    catch (const std::runtime_error &)
    {
        closedRejected = true;
    }
    require(closedRejected, "closed interface accepted a run");

    KernelInterface legacyClose;
    legacyClose.free();
    require(legacyClose.isClosed(), "legacy free did not forward to close");

    std::error_code ignored;
    std::filesystem::remove(malformed, ignored);
    std::filesystem::remove(unsupported, ignored);
    std::filesystem::remove(incomplete, ignored);
    return 0;
}
