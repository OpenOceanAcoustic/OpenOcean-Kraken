#ifndef OPEN_OCEAN_KRAKENC_INTERFACE_H
#define OPEN_OCEAN_KRAKENC_INTERFACE_H

#include "OpenOceanKrakencParams.h"
#include "OpenOceanKrakencSafety.h"

#include <memory>
#include <complex>
#include <string>
#include <vector>

namespace OpenOceanKrakenc
{
class ThreadPool;

class KernelInterface
{
public:
    KernelInterface();
    explicit KernelInterface(std::shared_ptr<ThreadPool> threadPool);
    [[deprecated("use KernelInterface(std::shared_ptr<ThreadPool>) or setProfileExecutor")]]
    explicit KernelInterface(ThreadPool &threadPool);
    ~KernelInterface();

    KernelInterface(const KernelInterface &) = delete;
    KernelInterface &operator=(const KernelInterface &) = delete;
    KernelInterface(KernelInterface &&) noexcept;
    KernelInterface &operator=(KernelInterface &&) noexcept;

    void setFieldThreads(int numThreads);
    int getFieldThreads() const;
    void setProfileExecutor(std::shared_ptr<ThreadPool> threadPool);
    void setThreadPool(std::shared_ptr<ThreadPool> threadPool);
    void detachThreadPool();

    [[deprecated("use setFieldThreads")]] void setNumThreads(int numThreads);
    [[deprecated("use getFieldThreads")]] int getNumThreads() const;
    int getHardwareThreads() const;
    [[deprecated("use setThreadPool(std::shared_ptr<ThreadPool>)")]]
    void setThreadPool(ThreadPool &threadPool);

    void runEigen();
    void runField();
    void run();

    void clearResults();
    void close();
    bool isClosed() const noexcept;
    [[deprecated("use close")]] void free();

    void set_Title(const std::string &title);
    void set_Freq(double freq);
    void set_freqvec(const Eigen::VectorXd &freqvec);
    void set_SSP(const std::vector<ssp::Range_Independent_Area> &sspInput);
    void set_AttenUnit(Atten_Mode mode);
    void set_Sz(const Eigen::VectorXd &values);
    void set_Sz(double start, double end, int count);
    void set_Rr(const Eigen::VectorXd &values);
    void set_Rr(double start, double end, int count);
    void set_Rz(const Eigen::VectorXd &values);
    void set_Rz(double start, double end, int count);
    void set_Ro(const Eigen::VectorXd &values);
    void set_Ro(double start, double end, int count);
    void set_cPhase(double cLow, double cHigh);
    void set_GridType(Grid_Mode type);
    void set_Rmax(double rMax);
    void set_SourceType(Source_Mode type);
    void set_MLimit(int limit);
    void set_CoherenceType(CoherenceType type);
    void set_ModeType(ModeType type);
    void set_ModeSampling(const Eigen::VectorXd &sourceDepths,
                          const Eigen::VectorXd &receiverDepths);
    void set_RunMode(Run_Mode mode);
    void set_Velocity_enable(bool enabled);
    void set_ReflCoef_Top(const std::vector<ReflectionCoef> &values);
    void set_ReflCoef_Bottom(const std::vector<ReflectionCoef> &values);
    void set_SBP(const Eigen::VectorXd &pattern,
                 const Eigen::VectorXd &anglesDegrees);

    ArrayView<const std::complex<float>> get_u_view(int sourceIndex) const;
    ArrayView<const std::complex<float>> get_v_view(int sourceIndex) const;
    ArrayView<const std::complex<float>> get_h_view(int sourceIndex) const;
    ArrayView<const std::complex<float>> get_u_view_all_sources() const;
    ArrayView<const std::complex<float>> get_v_view_all_sources() const;
    ArrayView<const std::complex<float>> get_h_view_all_sources() const;

    [[deprecated("use get_u_view")]] std::complex<float> *get_u(int sourceIndex);
    [[deprecated("use get_v_view")]] std::complex<float> *get_v(int sourceIndex);
    [[deprecated("use get_h_view")]] std::complex<float> *get_h(int sourceIndex);
    [[deprecated("use get_u_view_all_sources")]] std::complex<float> *get_u_AllSources();
    [[deprecated("use get_v_view_all_sources")]] std::complex<float> *get_v_AllSources();
    [[deprecated("use get_h_view_all_sources")]] std::complex<float> *get_h_AllSources();

    void export_result(const std::string &root);
    void export_mod(const std::string &root);
    void export_shd(const std::string &root, int dataType);

    LoadResult loadEnv(const std::string &envPath);
#if defined(OPENOCEAN_KRAKENC_LEGACY_JSON)
    LoadResult loadJson(const std::string &jsonPath);
#endif
    [[deprecated("use loadEnv to retain error details")]]
    bool from_env(const std::string &envPath);
#if defined(OPENOCEAN_KRAKENC_LEGACY_JSON)
    [[deprecated("use loadJson to retain error details")]]
    bool from_json(const std::string &jsonPath);
    bool to_json(const std::string &jsonPath) const;
    std::string to_json_string() const;
#endif

    OOKC_parameters &getParams();
    const OOKC_parameters &getParams() const;
    const OOKC_parameters &getParams_const() const;
    [[deprecated("use getOutput_const, getOutput_Copy, or result ArrayView APIs")]]
    OOKC_output &getOutput();
    const OOKC_output &getOutput() const;
    const OOKC_output &getOutput_const() const;
    OOKC_output getOutput_Copy() const;

private:
    class InterfaceImpl;

    void ensureAlive() const;
    void ensureResultsFresh() const;
    void invalidateResultViews() noexcept;

    std::unique_ptr<InterfaceImpl> impl_;
};
}

#endif
