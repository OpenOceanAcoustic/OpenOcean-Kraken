#include "OpenOceanKrakencKernelInterface.h"

#include <cassert>
#include <complex>
#include <limits>
#include <string>
#include <stdexcept>
#include <type_traits>
#include <vector>

using namespace OpenOceanKrakenc;

#undef assert
#define assert(condition) do { if (!(condition)) throw std::runtime_error("assertion failed: " #condition); } while (false)

namespace
{
template <typename T>
void instantiate()
{
    T value{};
    (void)value;
}
}

int main()
{
    instantiate<SSP_Mode>();
    instantiate<Media_Mode>();
    instantiate<AttenuationUnit>();
    instantiate<OceanAbsorptionModel>();
    instantiate<BC_Mode>();
    instantiate<Source_Mode>();
    instantiate<Grid_Mode>();
    instantiate<Run_Mode>();
    instantiate<CoherenceType>();
    instantiate<ModeType>();
    instantiate<Atten_Mode>();
    instantiate<HSInfo>();
    instantiate<Position>();
    instantiate<ReflectionCoef>();
    instantiate<InternalReflectionCoefInfo>();
    instantiate<ReflectionCoefInfo>();
    instantiate<SrcBmPat>();
    instantiate<BdryType>();
    instantiate<ssp::SSPLayer>();
    instantiate<ssp::SSPStructure>();
    instantiate<ssp::Range_Independent_Area>();
    instantiate<ssp::FlattenedData>();

    static_assert(std::is_same_v<decltype(&KernelInterface::from_env),
                                 bool (KernelInterface::*)(const std::string &)>);
    static_assert(std::is_same_v<decltype(&KernelInterface::from_json),
                                 bool (KernelInterface::*)(const std::string &)>);
    static_assert(std::is_same_v<decltype(&KernelInterface::to_json),
                                 bool (KernelInterface::*)(const std::string &) const>);
    static_assert(std::is_same_v<decltype(&KernelInterface::to_json_string),
                                 std::string (KernelInterface::*)() const>);
    static_assert(std::is_same_v<decltype(&KernelInterface::getHardwareThreads),
                                 int (KernelInterface::*)() const>);
    static_assert(std::is_same_v<decltype(&KernelInterface::getOutput_Copy),
                                 OOKC_output (KernelInterface::*)() const>);

    auto setTitle = &KernelInterface::set_Title;
    auto setFreq = &KernelInterface::set_Freq;
    auto setFreqVec = &KernelInterface::set_freqvec;
    auto setSsp = &KernelInterface::set_SSP;
    auto setAtten = &KernelInterface::set_AttenUnit;
    auto setPhase = &KernelInterface::set_cPhase;
    auto setGrid = &KernelInterface::set_GridType;
    auto setRange = &KernelInterface::set_Rmax;
    auto setSource = &KernelInterface::set_SourceType;
    auto setRun = &KernelInterface::set_RunMode;
    auto setVelocity = &KernelInterface::set_Velocity_enable;
    auto setTop = &KernelInterface::set_ReflCoef_Top;
    auto setBottom = &KernelInterface::set_ReflCoef_Bottom;
    auto setPattern = &KernelInterface::set_SBP;
    auto getU = &KernelInterface::get_u;
    auto getV = &KernelInterface::get_v;
    auto getH = &KernelInterface::get_h;
    auto getUAll = &KernelInterface::get_u_AllSources;
    auto getVAll = &KernelInterface::get_v_AllSources;
    auto getHAll = &KernelInterface::get_h_AllSources;
    auto exportResult = &KernelInterface::export_result;
    auto exportMod = &KernelInterface::export_mod;
    auto exportShd = &KernelInterface::export_shd;
    auto getParamsConst = &KernelInterface::getParams_const;
    auto getOutputConst = &KernelInterface::getOutput_const;
    (void)setTitle; (void)setFreq; (void)setFreqVec; (void)setSsp; (void)setAtten;
    (void)setPhase; (void)setGrid; (void)setRange; (void)setSource; (void)setRun;
    (void)setVelocity; (void)setTop; (void)setBottom; (void)setPattern;
    (void)getU; (void)getV; (void)getH; (void)getUAll; (void)getVAll; (void)getHAll;
    (void)exportResult; (void)exportMod; (void)exportShd;
    (void)getParamsConst; (void)getOutputConst;

    KernelInterface api;
    const Eigen::VectorXd values = Eigen::VectorXd::LinSpaced(3, 1.0, 3.0);
    api.set_Sz(values);
    api.set_Sz(10.0, 30.0, 3);
    api.set_Rr(values);
    api.set_Rr(100.0, 300.0, 3);
    api.set_Rz(values);
    api.set_Rz(20.0, 60.0, 3);
    api.set_Ro(values);
    api.set_Ro(0.0, 2.0, 3);
    api.set_Rmax(0.0);
    api.getParams().ReflectionCoef.isDeg = true;
    api.set_ReflCoef_Top({ReflectionCoef{30.0, 0.8, 180.0}});
    assert(!api.getParams_const().ReflectionCoef.isDeg);
    api.getParams().ReflectionCoef.isDeg = true;
    api.set_ReflCoef_Bottom({ReflectionCoef{30.0, 0.8, 180.0}});
    assert(!api.getParams_const().ReflectionCoef.isDeg);

    const OOKC_parameters &params = api.getParams_const();
    assert(params.Pos.NSz == 3);
    assert(params.Pos.NRr == 3);
    assert(params.Pos.NRz == 3);
    assert(params.Pos.NRo == 3);
    assert(params.Rmax == 0.0);
    assert(api.getHardwareThreads() >= 1);

    bool nonFinitePositionRejected = false;
    try
    {
        api.set_Rr(Eigen::VectorXd::Constant(
            1, std::numeric_limits<double>::infinity()));
    }
    catch (const std::invalid_argument &)
    {
        nonFinitePositionRejected = true;
    }
    assert(nonFinitePositionRejected);

    bool unorderedPatternRejected = false;
    try
    {
        api.set_SBP((Eigen::VectorXd(2) << 1.0, 1.0).finished(),
                    (Eigen::VectorXd(2) << 10.0, 5.0).finished());
    }
    catch (const std::invalid_argument &)
    {
        unorderedPatternRejected = true;
    }
    assert(unorderedPatternRejected);
    return 0;
}
