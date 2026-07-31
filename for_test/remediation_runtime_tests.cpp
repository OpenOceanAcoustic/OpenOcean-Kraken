#include "OpenOceanKrakencKernelInterface.h"
#include "ThreadPool.h"
#include "algorithm/AcousticCase.h"
#include "algorithm/ComplexMatrixBuilder.h"
#include "module/json_in_out.hpp"

#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <complex>
#include <filesystem>
#include <fstream>
#include <future>
#include <iostream>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using namespace OpenOceanKrakenc;
using namespace std::chrono_literals;

namespace
{
void require(bool condition, const std::string &message)
{
    if (!condition)
    {
        throw std::runtime_error(message);
    }
}

template <typename Exception, typename Function>
std::string captureException(Function function)
{
    try
    {
        function();
    }
    catch (const Exception &error)
    {
        return error.what();
    }
    throw std::runtime_error("expected exception was not thrown");
}

std::filesystem::path temporaryPath(
    const std::string &stem, const std::string &extension)
{
    static std::atomic<unsigned long> sequence{0};
    return std::filesystem::temp_directory_path() /
           (stem + "_" + std::to_string(sequence.fetch_add(1)) + extension);
}

AcousticCase validCase()
{
    AcousticCase input;
    input.title = "runtime-remediation";
    input.frequency = 50.0;
    input.referenceFrequency = 50.0;
    input.cLow = 1500.0;
    input.cHigh = 5000.0;
    input.rMaxKm = 0.0;
    input.top.type = AcousticBoundaryType::Vacuum;
    input.bottom.type = AcousticBoundaryType::Rigid;

    AcousticLayer layer;
    layer.baseMesh = 10;
    layer.topDepth = 0.0;
    layer.bottomDepth = 100.0;
    layer.samples = {
        {0.0, 1500.0, 0.0, 1.0, 0.0, 0.0},
        {100.0, 1500.0, 0.0, 1.0, 0.0, 0.0},
    };
    input.layers.push_back(layer);
    return input;
}

std::filesystem::path writeRmaxEnvironment(double rmax)
{
    const std::filesystem::path path =
        temporaryPath("ookc_alg004", ".env");
    std::ofstream stream(path);
    stream
        << "'ALG-004 runtime input'\n"
        << "50.0\n"
        << "1\n"
        << "'CVW'\n"
        << "10 0.0 100.0\n"
        << "0.0 1500.0 0.0 1.0 0.0 0.0\n"
        << "100.0 1500.0 /\n"
        << "'R' 0.0\n"
        << "1500.0 5000.0\n"
        << rmax << "\n"
        << "1\n"
        << "25.0 /\n"
        << "1\n"
        << "50.0 /\n";
    require(static_cast<bool>(stream), "failed to write ALG-004 ENV");
    return path;
}

struct RmaxObservation
{
    bool accepted = false;
    double value = 0.0;
};

RmaxObservation observeEnvRmax(double value)
{
    const std::filesystem::path path = writeRmaxEnvironment(value);
    std::error_code ignored;
    try
    {
        const AcousticCase parsed = readAcousticEnv(path);
        std::filesystem::remove(path, ignored);
        return {true, parsed.rMaxKm};
    }
    catch (const std::exception &)
    {
        std::filesystem::remove(path, ignored);
        return {};
    }
}

RmaxObservation observeJsonRmax(
    const OpenOcean_json &base, double value)
{
    OpenOcean_json document = base;
    document["Rmax"] = value;
    const std::filesystem::path path =
        temporaryPath("ookc_alg004", ".json");
    {
        std::ofstream stream(path);
        stream << document.dump(2) << '\n';
        require(static_cast<bool>(stream), "failed to write ALG-004 JSON");
    }
    OOKC_parameters parsed;
    const LoadResult loaded =
        read_json_file_result(path.string(), parsed);
    std::error_code ignored;
    std::filesystem::remove(path, ignored);
    return {loaded.ok, loaded.ok ? parsed.Rmax : 0.0};
}

RmaxObservation observeApiRmax(double value)
{
    KernelInterface api;
    try
    {
        api.set_Rmax(value);
        return {true, api.getParams_const().Rmax};
    }
    catch (const std::invalid_argument &)
    {
        return {};
    }
}

void testAlg004()
{
    const std::filesystem::path seedPath =
        std::filesystem::path(OPENOCEANKRAKENC_SOURCE_DIR) /
        "for_test/fixtures/two_profile_small.env";
    KernelInterface seed;
    const LoadResult seedLoaded = seed.loadEnv(seedPath.string());
    require(seedLoaded.ok,
            "ALG-004 checked-in ENV/FLP seed must load");
    const OpenOcean_json base =
        OpenOcean_json::parse(seed.to_json_string());

    for (const double input : {-1.0, 0.0, 10.0})
    {
        const bool expectedAccepted = input >= 0.0;
        const std::array<std::pair<const char *, RmaxObservation>, 3>
            observations = {{
                {"ENV", observeEnvRmax(input)},
                {"JSON", observeJsonRmax(base, input)},
                {"API", observeApiRmax(input)},
            }};
        for (const auto &route : observations)
        {
            require(
                route.second.accepted == expectedAccepted,
                std::string(route.first) +
                    " RMax acceptance differs from the unified contract");
            if (expectedAccepted)
            {
                require(
                    route.second.value == input,
                    std::string(route.first) +
                        " changed an accepted RMax value");
            }
        }
    }

    const double infinity =
        std::numeric_limits<double>::infinity();
    require(!observeEnvRmax(infinity).accepted,
            "ENV must reject non-finite RMax");
    require(!observeApiRmax(infinity).accepted,
            "API must reject non-finite RMax");
}

void expectMatrixRejected(
    AcousticCase input, const std::string &fragment)
{
    const std::string message =
        captureException<std::invalid_argument>([&] {
            static_cast<void>(buildAcousticMatrix(input, 1));
        });
    require(message.find(fragment) != std::string::npos,
            "ALG-007 error did not contain: " + fragment);
}

void testAlg007()
{
    {
        AcousticCase input = validCase();
        input.layers.front().samples.front().depth = 1.0;
        expectMatrixRejected(input, "span");
    }
    {
        AcousticCase input = validCase();
        input.layers.front().samples.clear();
        expectMatrixRejected(input, "at least two");
    }
    {
        AcousticCase input = validCase();
        input.layers.front().samples.resize(1);
        expectMatrixRejected(input, "at least two");
    }
    {
        AcousticCase input = validCase();
        AcousticSample interior = input.layers.front().samples.front();
        input.layers.front().samples.insert(
            input.layers.front().samples.end() - 1, interior);
        expectMatrixRejected(input, "strictly increasing");
    }
    {
        AcousticCase input = validCase();
        AcousticSample interior = input.layers.front().samples.front();
        interior.depth = -1.0;
        input.layers.front().samples.insert(
            input.layers.front().samples.end() - 1, interior);
        expectMatrixRejected(input, "strictly increasing");
    }

    const std::array<double AcousticSample::*, 6> finiteMembers = {
        &AcousticSample::depth,
        &AcousticSample::cp,
        &AcousticSample::cs,
        &AcousticSample::rho,
        &AcousticSample::alphaP,
        &AcousticSample::alphaS,
    };
    for (double AcousticSample::*member : finiteMembers)
    {
        AcousticCase input = validCase();
        std::size_t sampleIndex = 0;
        if (member == &AcousticSample::depth)
        {
            AcousticSample interior = input.layers.front().samples.front();
            interior.depth = 50.0;
            input.layers.front().samples.insert(
                input.layers.front().samples.end() - 1, interior);
            sampleIndex = 1;
        }
        input.layers.front().samples[sampleIndex].*member =
            std::numeric_limits<double>::quiet_NaN();
        expectMatrixRejected(input, "finite");
    }

    {
        AcousticCase input = validCase();
        input.layers.front().topDepth =
            std::numeric_limits<double>::quiet_NaN();
        expectMatrixRejected(input, "bounds");
    }
    {
        AcousticCase input = validCase();
        input.layers.front().bottomDepth =
            std::numeric_limits<double>::infinity();
        expectMatrixRejected(input, "bounds");
    }
    {
        AcousticCase input = validCase();
        input.layers.front().bottomDepth =
            input.layers.front().topDepth;
        expectMatrixRejected(input, "greater than top");
    }
    {
        AcousticCase input = validCase();
        input.layers.front().bottomDepth = -1.0;
        expectMatrixRejected(input, "greater than top");
    }
    {
        AcousticCase input = validCase();
        input.layers.front().samples.back().depth = 99.0;
        expectMatrixRejected(input, "span");
    }

    const AcousticMatrix matrix =
        buildAcousticMatrix(validCase(), 1);
    require(!matrix.cp.empty(),
            "valid matrix construction regressed");
}

void testAlg019Completion()
{
    ThreadPool pool(1);
    auto task = pool.enqueue([&pool] {
        const std::string message =
            captureException<std::logic_error>([&] {
                pool.wait_completion();
            });
        return message.find("current task") != std::string::npos;
    });
    require(task.get(),
            "pool-internal wait_completion must reject the current task");
    require(pool.wait_completion_for(100ms),
            "external completion wait must still succeed");
}

void testAlg019Id()
{
    {
        ThreadPool pool(1);
        auto task = pool.enqueue_with_id("self", [&pool] {
            const std::string message =
                captureException<std::logic_error>([&] {
                    pool.wait_id("self");
                });
            return message.find("current task ID") !=
                   std::string::npos;
        });
        require(task.get(),
                "identified task must reject waiting for its own ID");
        require(pool.wait_id_for("self", 100ms),
                "external ID wait must still succeed");
    }

    {
        ThreadPool pool(2);
        std::promise<void> otherStartedPromise;
        std::shared_future<void> otherStarted =
            otherStartedPromise.get_future().share();
        std::promise<void> releaseOtherPromise;
        std::shared_future<void> releaseOther =
            releaseOtherPromise.get_future().share();

        auto other = pool.enqueue_with_id(
            "other", [&] {
                otherStartedPromise.set_value();
                releaseOther.wait();
            });
        auto waiter = pool.enqueue_with_id(
            "waiter", [&] {
                otherStarted.wait();
                pool.wait_id("other");
                return true;
            });

        otherStarted.wait();
        releaseOtherPromise.set_value();
        other.get();
        require(waiter.get(),
                "worker must be allowed to wait for a different ID");
    }
}

void testAlg019Timed()
{
    ThreadPool pool(1);
    auto completion = pool.enqueue([&pool] {
        const auto started =
            std::chrono::steady_clock::now();
        const bool result =
            pool.wait_completion_for(250ms);
        const auto elapsed =
            std::chrono::steady_clock::now() - started;
        return !result && elapsed < 100ms;
    });
    require(completion.get(),
            "timed completion self-wait must return false within 100 ms");

    auto identified = pool.enqueue_with_id(
        "timed-self", [&pool] {
            const auto started =
                std::chrono::steady_clock::now();
            const bool result =
                pool.wait_id_for("timed-self", 250ms);
            const auto elapsed =
                std::chrono::steady_clock::now() - started;
            return !result && elapsed < 100ms;
        });
    require(identified.get(),
            "timed ID self-wait must return false within 100 ms");
}

using ScenarioFunction = void (*)();

std::map<std::string, ScenarioFunction> &scenarioRegistry()
{
    static std::map<std::string, ScenarioFunction> scenarios;
    return scenarios;
}

struct RegisterScenario
{
    RegisterScenario(std::string name, ScenarioFunction function)
    {
        const bool inserted =
            scenarioRegistry().emplace(std::move(name), function).second;
        require(inserted, "duplicate remediation runtime scenario");
    }
};

const RegisterScenario registerAlg004{"ALG-004", &testAlg004};
const RegisterScenario registerAlg007{"ALG-007", &testAlg007};
const RegisterScenario registerAlg019Completion{
    "ALG-019-completion", &testAlg019Completion};
const RegisterScenario registerAlg019Id{
    "ALG-019-id", &testAlg019Id};
const RegisterScenario registerAlg019Timed{
    "ALG-019-timed", &testAlg019Timed};
}

int main(int argc, char **argv)
{
    if (argc != 2)
    {
        std::cerr << "exactly one remediation scenario is required\n";
        return 2;
    }
    const auto found = scenarioRegistry().find(argv[1]);
    if (found == scenarioRegistry().end())
    {
        std::cerr << "unknown remediation scenario: " << argv[1] << '\n';
        return 2;
    }
    try
    {
        found->second();
        return 0;
    }
    catch (const std::exception &error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
