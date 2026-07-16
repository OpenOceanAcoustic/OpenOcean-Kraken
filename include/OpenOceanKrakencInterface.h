#ifndef OPEN_OCEAN_KRAKENC_INTERFACE_H
#define OPEN_OCEAN_KRAKENC_INTERFACE_H

#include "OpenOceanKrakencParams.h"

#include <memory>
#include <complex>
#include <string>
#include <vector>

class ThreadPool;

namespace OpenOceanKrakenc
{
class Interface
{
public:
    Interface();
    explicit Interface(ThreadPool &threadPool);
    ~Interface();

    Interface(const Interface &) = delete;
    Interface &operator=(const Interface &) = delete;
    Interface(Interface &&) noexcept;
    Interface &operator=(Interface &&) noexcept;

    void setNumThreads(int numThreads);
    int getNumThreads() const;
    int getHardwareThreads() const;
    void setThreadPool(ThreadPool &threadPool);

    void runEigen();
    void runField();
    void run();

    void clearResults();
    void free();

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
    void set_RunMode(Run_Mode mode);
    void set_Velocity_enable(bool enabled);
    void set_ReflCoef_Top(const std::vector<ReflectionCoef> &values);
    void set_ReflCoef_Bottom(const std::vector<ReflectionCoef> &values);
    void set_SBP(const Eigen::VectorXd &pattern,
                 const Eigen::VectorXd &anglesDegrees);

    std::complex<float> *get_u(int sourceIndex);
    std::complex<float> *get_v(int sourceIndex);
    std::complex<float> *get_h(int sourceIndex);
    std::complex<float> *get_u_AllSources();
    std::complex<float> *get_v_AllSources();
    std::complex<float> *get_h_AllSources();

    void export_result(const std::string &root);
    void export_mod(const std::string &root);
    void export_shd(const std::string &root, int dataType);

    bool from_env(const std::string &envPath);
    bool from_json(const std::string &jsonPath);
    bool to_json(const std::string &jsonPath) const;
    std::string to_json_string() const;

    OOKC_parameters &getParams();
    const OOKC_parameters &getParams() const;
    const OOKC_parameters &getParams_const() const;
    OOKC_output &getOutput();
    const OOKC_output &getOutput() const;
    const OOKC_output &getOutput_const() const;
    OOKC_output getOutput_Copy() const;

private:
    class InterfaceImpl;

    void ensureAlive() const;
    void ensureResultsFresh() const;

    std::unique_ptr<InterfaceImpl> impl_;
};
}

#endif
