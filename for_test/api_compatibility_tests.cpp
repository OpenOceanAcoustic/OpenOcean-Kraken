#include "OpenOceanKrakencInterface.h"

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

    static_assert(std::is_same_v<decltype(&Interface::from_env),
                                 bool (Interface::*)(const std::string &)>);
    static_assert(std::is_same_v<decltype(&Interface::from_json),
                                 bool (Interface::*)(const std::string &)>);
    static_assert(std::is_same_v<decltype(&Interface::to_json),
                                 bool (Interface::*)(const std::string &) const>);
    static_assert(std::is_same_v<decltype(&Interface::to_json_string),
                                 std::string (Interface::*)() const>);
    static_assert(std::is_same_v<decltype(&Interface::getHardwareThreads),
                                 int (Interface::*)() const>);
    static_assert(std::is_same_v<decltype(&Interface::getOutput_Copy),
                                 OOKC_output (Interface::*)() const>);

    auto setTitle = &Interface::set_Title;
    auto setFreq = &Interface::set_Freq;
    auto setFreqVec = &Interface::set_freqvec;
    auto setSsp = &Interface::set_SSP;
    auto setAtten = &Interface::set_AttenUnit;
    auto setPhase = &Interface::set_cPhase;
    auto setGrid = &Interface::set_GridType;
    auto setRange = &Interface::set_Rmax;
    auto setSource = &Interface::set_SourceType;
    auto setRun = &Interface::set_RunMode;
    auto setVelocity = &Interface::set_Velocity_enable;
    auto setTop = &Interface::set_ReflCoef_Top;
    auto setBottom = &Interface::set_ReflCoef_Bottom;
    auto setPattern = &Interface::set_SBP;
    auto getU = &Interface::get_u;
    auto getV = &Interface::get_v;
    auto getH = &Interface::get_h;
    auto getUAll = &Interface::get_u_AllSources;
    auto getVAll = &Interface::get_v_AllSources;
    auto getHAll = &Interface::get_h_AllSources;
    auto exportResult = &Interface::export_result;
    auto exportMod = &Interface::export_mod;
    auto exportShd = &Interface::export_shd;
    auto getParamsConst = &Interface::getParams_const;
    auto getOutputConst = &Interface::getOutput_const;
    (void)setTitle; (void)setFreq; (void)setFreqVec; (void)setSsp; (void)setAtten;
    (void)setPhase; (void)setGrid; (void)setRange; (void)setSource; (void)setRun;
    (void)setVelocity; (void)setTop; (void)setBottom; (void)setPattern;
    (void)getU; (void)getV; (void)getH; (void)getUAll; (void)getVAll; (void)getHAll;
    (void)exportResult; (void)exportMod; (void)exportShd;
    (void)getParamsConst; (void)getOutputConst;

    Interface api;
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
