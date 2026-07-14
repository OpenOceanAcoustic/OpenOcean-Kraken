#include "test.h"
#include "AttenMod.h"
#include "BCImpedanceMod.h"
#include "OpenOceanKrakenInterface.h"
#include "pchipMod.h"
#include "RefCoef.h"
#include "Scatter.h"
#include "splinec.h"
#include "sspMod.h"
#include "ThreadPool.h"
#include "util.h"
#include "field.h"
#include "input_Freq.hpp"
#include "input_reflcoef.hpp"
#include "input_sbp.hpp"
#include "input_SSP.hpp"
#include "input_Sz_Rz_RR.hpp"

#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <vector>

namespace OpenOceanKraken
{
    bool read_refCoef_file(const std::string &envPath, OOK_parameters &params, std::string pattern);
}

namespace fs = std::filesystem;
using OpenOceanKraken::Atten_Mode;
using OpenOceanKraken::AttenuationUnit;
using OpenOceanKraken::BCImpedance;
using OpenOceanKraken::BC_Mode;
using OpenOceanKraken::CoherenceType;
using OpenOceanKraken::CRCI;
using OpenOceanKraken::CSpline;
using OpenOceanKraken::ElasticDN;
using OpenOceanKraken::ElasticUP;
using OpenOceanKraken::EigenParams;
using OpenOceanKraken::Evaluate;
using OpenOceanKraken::Grid_Mode;
using OpenOceanKraken::HSInfo;
using OpenOceanKraken::Interface;
using OpenOceanKraken::InterpolateIRC;
using OpenOceanKraken::InterpolateReflectionCoefficient;
using OpenOceanKraken::InternalReflectionCoefInfo;
using OpenOceanKraken::KupIng;
using OpenOceanKraken::Media_Mode;
using OpenOceanKraken::ModeType;
using OpenOceanKraken::OceanAbsorptionModel;
using OpenOceanKraken::OOK_output;
using OpenOceanKraken::OOK_parameters;
using OpenOceanKraken::PCHIP;
using OpenOceanKraken::ReflectionCoef;
using OpenOceanKraken::Run_Mode;
using OpenOceanKraken::ScatterRoot;
using OpenOceanKraken::SplineALL;
using OpenOceanKraken::Source_Mode;
using OpenOceanKraken::SSP_Mode;
using OpenOceanKraken::TridMtx;
using OpenOceanKraken::VSpline;
using OpenOceanKraken::addOceanAbsorption;
using OpenOceanKraken::field;
using OpenOceanKraken::fprime_interior;
using OpenOceanKraken::fprime_left_end;
using OpenOceanKraken::fprime_right_end;
using OpenOceanKraken::getSourceEnv;
using OpenOceanKraken::h_del;
using OpenOceanKraken::input_Freq;
using OpenOceanKraken::input_reflcoef;
using OpenOceanKraken::input_sbp;
using OpenOceanKraken::input_SSP;
using OpenOceanKraken::input_Sz_Rz_RR;
using OpenOceanKraken::outputBase;
using OpenOceanKraken::paramsBase;
using OpenOceanKraken::parseAttenuation;
using OpenOceanKraken::read_refCoef_file;
using OpenOceanKraken::spline;
using OpenOceanKraken::ssp::Range_Independent_Area;
namespace Util = OpenOceanKraken::Util;

namespace
{
    std::vector<std::string> expected_env_failures()
    {
        return {"neggradC_brc.env", "neggradK_brc.env", "neggradC_irc.env", "neggradK_irc.env"};
    }

    bool is_expected_env_failure(const std::string &name)
    {
        const auto failures = expected_env_failures();
        return std::find(failures.begin(), failures.end(), name) != failures.end();
    }

    fs::path positive_fixture_root(const fs::path &test_root)
    {
        const fs::path good_examples = test_root / "good_examples";
        return fs::exists(good_examples) ? good_examples : test_root;
    }

    void require_position_consistency(ook_test::TestRunner &test, const OOK_parameters &params, const std::string &case_name)
    {
        test.require(params.Pos.NSz > 0, case_name + ": NSz must be positive");
        test.require(params.Pos.NRz > 0, case_name + ": NRz must be positive");
        test.require(params.Pos.NRr > 0, case_name + ": NRr must be positive");
        test.require(params.Pos.NRo > 0, case_name + ": NRo must be positive");
        test.require(params.Pos.Sz.size() == params.Pos.NSz, case_name + ": Sz size mismatch");
        test.require(params.Pos.Rz.size() == params.Pos.NRz, case_name + ": Rz size mismatch");
        test.require(params.Pos.Rr.size() == params.Pos.NRr, case_name + ": Rr size mismatch");
        test.require(params.Pos.Ro.size() == params.Pos.NRo, case_name + ": Ro size mismatch");
        test.require(params.Pos.GridType == Grid_Mode::MODE_R_Rectangular, case_name + ": FLP should select rectangular grid");
    }

    void require_ssp_consistency(ook_test::TestRunner &test, const OOK_parameters &params, const std::string &case_name)
    {
        test.require(!params.sspInput.empty(), case_name + ": sspInput must not be empty");
        test.require(params.SSP.size() == params.sspInput.size(), case_name + ": SSP conversion count mismatch");
        for (size_t i = 0; i < params.sspInput.size(); ++i)
        {
            const auto &input_area = params.sspInput[i];
            const auto &ssp = params.SSP[i];
            test.require(!input_area.layers.empty(), case_name + ": input area has no layers");
            test.require(ssp.NMedia == static_cast<int>(input_area.layers.size()), case_name + ": NMedia mismatch");
            test.require(ssp.NPts.size() == ssp.NMedia, case_name + ": NPts size mismatch");
            test.require(ssp.NMesh.size() == ssp.NMedia, case_name + ": NMesh size mismatch");
            test.require(ssp.alphaR.size() == ssp.rho.size(), case_name + ": alphaR/rho size mismatch");
        }
    }

    void test_lifecycle_and_guards(ook_test::TestRunner &test)
    {
        Interface iface;
        test.require(iface.getNumThreads() == 1, "Default Interface must use one thread");
        test.requireThrows([&]() { iface.setNumThreads(0); }, "setNumThreads(0) must throw");
        test.requireThrows([&]() { iface.runEigen(); }, "runEigen without ThreadPool must throw");
        test.requireThrows([&]() { iface.getOutput_Copy(); }, "getOutput_Copy must be disabled");

        Interface freed;
        freed.free();
        freed.free();
        test.requireThrows([&]() { freed.getParams(); }, "getParams after free must throw");
    }

    void test_top_bottom_line(ook_test::TestRunner &test)
    {
        Range_Independent_Area area;
        area.set_Bottom_Line(5000.0, 1600.0, 0.8, 0.0, 0.0, 1.8);
        area.set_Top_Line(0.0, 1500.0, 0.0, 0.0, 0.0, 1.0);
        test.requireNear(area.HSTop.Depth, 0.0, 1.0e-9, "Top depth must be written to HSTop");
        test.requireNear(area.HSTop.alphaR, 1500.0, 1.0e-9, "Top alphaR must be written to HSTop");
        test.requireNear(area.HSBot.Depth, 5000.0, 1.0e-9, "Top setter must not overwrite HSBot depth");
        test.requireNear(area.HSBot.alphaR, 1600.0, 1.0e-9, "Top setter must not overwrite HSBot alphaR");
        test.require(area.HSTop.BC == BC_Mode::MODE_A_Half_space, "Top setter must set top BC");
        test.require(area.HSBot.BC == BC_Mode::MODE_A_Half_space, "Bottom setter must keep bottom BC");
    }

    Eigen::MatrixXcd spline_fixture(int n)
    {
        Eigen::MatrixXcd c(4, n);
        c.setZero();
        for (int i = 0; i < n; ++i)
        {
            c(0, i) = std::complex<double>(1.0 + i, 0.25 * i);
            c(1, i) = std::complex<double>(0.5 + 0.1 * i, 0.05 * i);
        }
        return c;
    }

    void test_interpolation_helpers(ook_test::TestRunner &test)
    {
        Eigen::VectorXd x2(2);
        x2 << 0.0, 1.0;
        Eigen::VectorXcd y2(2);
        y2 << std::complex<double>(1.0, 0.0), std::complex<double>(3.0, 2.0);
        Eigen::MatrixXcd poly2(4, 2);
        Eigen::MatrixXcd work2(4, 2);
        poly2.setZero();
        work2.setZero();
        PCHIP(x2, y2, 2, poly2, work2);
        test.requireNear(std::real(poly2(1, 0)), 2.0, 1.0e-12, "PCHIP N=2 slope real mismatch");
        test.requireNear(std::imag(poly2(1, 0)), 2.0, 1.0e-12, "PCHIP N=2 slope imag mismatch");

        Eigen::VectorXd x4(4);
        x4 << 0.0, 1.0, 2.0, 4.0;
        Eigen::VectorXcd y4(4);
        y4 << std::complex<double>(0.0, 0.0),
            std::complex<double>(1.0, -1.0),
            std::complex<double>(1.5, -0.5),
            std::complex<double>(4.0, -2.0);
        Eigen::MatrixXcd poly4(4, 4);
        Eigen::MatrixXcd work4(4, 4);
        poly4.setZero();
        work4.setZero();
        PCHIP(x4, y4, 4, poly4, work4);
        test.require(std::abs(poly4(0, 0)) < 1.0e-12, "PCHIP must retain first sample");

        test.requireNear(fprime_interior(1.0, 2.0, 10.0), 3.0, 1.0e-12, "positive interior derivative clamp");
        test.requireNear(fprime_interior(-1.0, -2.0, -10.0), -3.0, 1.0e-12, "negative interior derivative clamp");
        test.requireNear(fprime_interior(1.0, -2.0, 10.0), 0.0, 1.0e-12, "mixed-sign interior derivative clamp");
        test.requireNear(fprime_left_end(1.0, 2.0, -1.0), 0.0, 1.0e-12, "left-end sign guard");
        test.requireNear(fprime_left_end(1.0, -2.0, 10.0), 3.0, 1.0e-12, "left-end monotonic clamp");
        test.requireNear(fprime_right_end(1.0, 2.0, -1.0), 0.0, 1.0e-12, "right-end sign guard");
        test.requireNear(fprime_right_end(-1.0, 2.0, 10.0), 6.0, 1.0e-12, "right-end monotonic clamp");

        double h1 = 0.0;
        double h2 = 0.0;
        std::complex<double> del1;
        std::complex<double> del2;
        h_del(x4, y4, 1, h1, h2, del1, del2);
        test.requireNear(h1, 1.0, 1.0e-12, "h_del h1 mismatch");
        test.requireNear(h2, 1.0, 1.0e-12, "h_del h2 mismatch");

        auto c2 = spline_fixture(2);
        auto c3 = spline_fixture(3);
        auto c4 = spline_fixture(4);
        CSpline(x2, c2, 2, 0, 0, 2);
        CSpline(x4, c4, 4, 0, 0, 4);
        CSpline(x4, c4, 4, 1, 1, 4);
        CSpline(x4, c4, 4, 2, 2, 4);
        Eigen::VectorXd x3(3);
        x3 << 0.0, 1.0, 2.0;
        CSpline(x3, c3, 3, 0, 0, 3);

        Eigen::VectorXcd coeffs(8);
        coeffs << std::complex<double>(1.0, 0.0), std::complex<double>(2.0, 0.0),
            std::complex<double>(3.0, 0.0), std::complex<double>(4.0, 0.0),
            std::complex<double>(2.0, 0.0), std::complex<double>(1.0, 0.0),
            std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0);
        Eigen::VectorXcd eval_points(2);
        eval_points << std::complex<double>(0.25, 0.0), std::complex<double>(0.75, 0.0);
        VSpline(x2, coeffs, 2, 2, eval_points, 2);
        test.require(std::real(eval_points(0)) > 1.0, "VSpline must evaluate first interval");

        const std::complex<double> raw_coeffs[4] = {
            std::complex<double>(1.0, 0.0),
            std::complex<double>(2.0, 0.0),
            std::complex<double>(3.0, 0.0),
            std::complex<double>(4.0, 0.0)};
        test.requireNear(std::real(spline(raw_coeffs, 0.5)), 2.458333333333333, 1.0e-12, "spline polynomial mismatch");

        int seg = 0;
        double dist = 0.25;
        std::complex<double> f;
        std::complex<double> fx;
        std::complex<double> fxx;
        SplineALL(c4, seg, dist, f, fx, fxx);
        test.require(std::abs(f) > 0.0, "SplineALL must produce a value");
    }

    void test_attenuation_and_reflection_helpers(ook_test::TestRunner &test)
    {
        auto unit = [](AttenuationUnit attn, OceanAbsorptionModel model = OceanAbsorptionModel::None) {
            Atten_Mode mode;
            mode.attnUnit = attn;
            mode.absModel = model;
            return mode;
        };

        test.requireNear(parseAttenuation(100.0, 100.0, 0.2, 200.0, 2.0, 1500.0, unit(AttenuationUnit::MODE_N_Nepers_per_m)), 0.2, 1.0e-12, "N attenuation mismatch");
        test.require(parseAttenuation(100.0, 100.0, 0.2, 200.0, 2.0, 1500.0, unit(AttenuationUnit::MODE_M_dB_per_m)) > 0.0, "M attenuation must be positive");
        test.require(parseAttenuation(100.0, 100.0, 0.2, 200.0, 2.0, 1500.0, unit(AttenuationUnit::MODE_m_dB_per_m)) > 0.0, "m attenuation low-frequency branch");
        test.require(parseAttenuation(300.0, 100.0, 0.2, 200.0, 2.0, 1500.0, unit(AttenuationUnit::MODE_m_dB_per_m)) > 0.0, "m attenuation high-frequency branch");
        test.require(parseAttenuation(100.0, 100.0, 0.2, 200.0, 2.0, 1500.0, unit(AttenuationUnit::MODE_F_dB_per_m_kHz)) > 0.0, "F attenuation must be positive");
        test.requireNear(parseAttenuation(100.0, 100.0, 0.2, 200.0, 2.0, 0.0, unit(AttenuationUnit::MODE_W_db_per_lambda)), 0.0, 1.0e-12, "W attenuation zero-c guard");
        test.require(parseAttenuation(100.0, 100.0, 0.2, 200.0, 2.0, 1500.0, unit(AttenuationUnit::MODE_Q_Quality_Factor)) > 0.0, "Q attenuation must be positive");
        test.requireNear(parseAttenuation(100.0, 100.0, 0.0, 200.0, 2.0, 1500.0, unit(AttenuationUnit::MODE_Q_Quality_Factor)), 0.0, 1.0e-12, "Q attenuation zero-alpha guard");
        test.require(parseAttenuation(100.0, 100.0, 0.2, 200.0, 2.0, 1500.0, unit(AttenuationUnit::MODE_L_params_lose)) > 0.0, "L attenuation must be positive");
        test.require(addOceanAbsorption(0.0, 1000.0, unit(AttenuationUnit::MODE_W_db_per_lambda, OceanAbsorptionModel::Thorpe)) > 0.0, "Thorpe absorption must be positive");
        test.require(addOceanAbsorption(0.0, 1000.0, unit(AttenuationUnit::MODE_W_db_per_lambda, OceanAbsorptionModel::FrancGarr)) > 0.0, "Franc-Garr absorption must be positive");

        double z = 0.0;
        double c = 1500.0;
        double alpha = 0.1;
        double freq = 100.0;
        double freq0 = 100.0;
        double beta = 1.0;
        double ft = 200.0;
        auto crci = CRCI(z, c, alpha, freq, freq0, unit(AttenuationUnit::MODE_W_db_per_lambda), beta, ft);
        test.requireNear(std::real(crci), 1500.0, 1.0e-12, "CRCI real component mismatch");

        Eigen::Matrix<ReflectionCoef, 1, Eigen::Dynamic> table(3);
        table(0).theta = 0.0;
        table(0).R = 0.2;
        table(0).phi = 1.0;
        table(1).theta = 10.0;
        table(1).R = 0.4;
        table(1).phi = 2.0;
        table(2).theta = 20.0;
        table(2).R = 0.6;
        table(2).phi = 3.0;
        ReflectionCoef interp{-1.0, 1.0, 1.0};
        InterpolateReflectionCoefficient(interp, table);
        test.requireNear(interp.R, 0.0, 1.0e-12, "reflection lower bound must zero");
        interp.theta = 25.0;
        interp.R = 1.0;
        InterpolateReflectionCoefficient(interp, table);
        test.requireNear(interp.R, 0.0, 1.0e-12, "reflection upper bound must zero");
        interp.theta = 15.0;
        InterpolateReflectionCoefficient(interp, table);
        test.requireNear(interp.R, 0.5, 1.0e-12, "reflection interpolation mismatch");

        InternalReflectionCoefInfo bad_irc;
        std::complex<double> f_out;
        std::complex<double> g_out;
        int power = 0;
        test.requireThrows([&]() { InterpolateIRC({1.0, 0.0}, f_out, g_out, power, bad_irc); }, "unset IRC must throw");

        InternalReflectionCoefInfo one;
        one.isSet = true;
        one.xTab.resize(1);
        one.fTab.resize(1);
        one.gTab.resize(1);
        one.iTab.resize(1);
        one.xTab << 10.0;
        one.fTab << std::complex<double>(2.0, 1.0);
        one.gTab << std::complex<double>(3.0, -1.0);
        one.iTab << 4;
        InterpolateIRC({5.0, 0.0}, f_out, g_out, power, one);
        test.require(power == 4, "single-point IRC power mismatch");

        InternalReflectionCoefInfo irc;
        irc.isSet = true;
        irc.xTab.resize(3);
        irc.fTab.resize(3);
        irc.gTab.resize(3);
        irc.iTab.resize(3);
        irc.xTab << 100.0, 200.0, 300.0;
        irc.fTab << std::complex<double>(1.0, 0.0), std::complex<double>(2.0, 0.0), std::complex<double>(3.0, 0.0);
        irc.gTab << std::complex<double>(2.0, 0.0), std::complex<double>(4.0, 0.0), std::complex<double>(6.0, 0.0);
        irc.iTab << 0, 1, 2;
        InterpolateIRC({250.0, 0.0}, f_out, g_out, power, irc);
        test.require(power == 1, "multi-point IRC interpolation power mismatch");
    }

    OOK_parameters minimal_boundary_params(BC_Mode top_bc, BC_Mode bot_bc)
    {
        OOK_parameters params;
        params.freqinfo.freq = 50.0;
        params.SSP.resize(1);
        auto &ssp = params.SSP[0];
        ssp.NMedia = 1;
        ssp.FirstAcoustic = 0;
        ssp.LastAcoustic = 0;
        ssp.NPts.resize(1);
        ssp.NMesh.resize(1);
        ssp.offset.resize(1);
        ssp.interp_offset.resize(1);
        ssp.depth.resize(1);
        ssp.NPts << 2;
        ssp.NMesh << 1;
        ssp.offset << 0;
        ssp.interp_offset << 0;
        ssp.depth << 100.0;
        ssp.HSTop.BC = top_bc;
        ssp.HSTop.cp = {1500.0, 0.0};
        ssp.HSTop.cs = {0.0, 0.0};
        ssp.HSTop.rho = 1.0;
        ssp.HSBot.BC = bot_bc;
        ssp.HSBot.cp = {1600.0, 0.0};
        ssp.HSBot.cs = {0.0, 0.0};
        ssp.HSBot.rho = 1.5;
        return params;
    }

    TridMtx minimal_trid()
    {
        TridMtx trid;
        trid.resize(4, 1, 1);
        trid.N << 2;
        trid.Loc << 0;
        trid.h << 1.0;
        trid.B1.setConstant(0.1);
        trid.B2.setConstant(0.2);
        trid.B3.setConstant(0.3);
        trid.B4.setConstant(0.4);
        trid.rho.setConstant(1.2);
        return trid;
    }

    void test_boundary_and_scatter_helpers(ook_test::TestRunner &test)
    {
        test.require(std::real(ScatterRoot({4.0, 0.0})) > 0.0, "ScatterRoot positive branch");
        test.require(std::imag(ScatterRoot({-4.0, 0.0})) < 0.0, "ScatterRoot negative branch");
        test.requireNear(std::abs(KupIng(0.0, {1.0, 0.0}, 1.0, {2.0, 0.0}, 2.0, {1.0, 0.0}, {1.0, 0.0})), 0.0, 1.0e-12, "KupIng zero sigma branch");
        test.require(std::abs(KupIng(0.5, {1.0, 0.0}, 1.0, {2.0, 0.0}, 2.0, {1.0, 0.0}, {0.5, 0.0})) > 0.0, "KupIng roughness branch");

        std::complex<double> f;
        std::complex<double> g;
        int iPower = 0;
        int modeCount = 0;

        auto vacuum_params = minimal_boundary_params(BC_Mode::MODE_V_Vacuum, BC_Mode::MODE_R_Rigid);
        auto trid = minimal_trid();
        BCImpedance(0, 0.01, true, f, g, iPower, true, trid, vacuum_params, modeCount);
        test.requireNear(std::real(f), 1.0, 1.0e-12, "vacuum top f mismatch");
        BCImpedance(0, 0.01, false, f, g, iPower, true, trid, vacuum_params, modeCount);
        test.requireNear(std::real(g), 1.0, 1.0e-12, "rigid bottom g mismatch");

        auto half_params = minimal_boundary_params(BC_Mode::MODE_A_Half_space, BC_Mode::MODE_A_Half_space);
        trid = minimal_trid();
        BCImpedance(0, 0.01, true, f, g, iPower, false, trid, half_params, modeCount);
        test.require(std::isfinite(std::real(f)), "half-space acoustic top must be finite");
        half_params.SSP[0].HSBot.cs = {800.0, 0.0};
        trid = minimal_trid();
        BCImpedance(0, 0.01, false, f, g, iPower, true, trid, half_params, modeCount);
        test.require(std::isfinite(std::real(g)), "half-space elastic bottom must be finite");

        auto file_params = minimal_boundary_params(BC_Mode::MODE_F_File, BC_Mode::MODE_F_File);
        file_params.ReflectionCoef.RTop.resize(2);
        file_params.ReflectionCoef.RBot.resize(2);
        for (int i = 0; i < 2; ++i)
        {
            file_params.ReflectionCoef.RTop(i).theta = i * 90.0;
            file_params.ReflectionCoef.RTop(i).R = 0.1 + 0.1 * i;
            file_params.ReflectionCoef.RTop(i).phi = 0.0;
            file_params.ReflectionCoef.RBot(i) = file_params.ReflectionCoef.RTop(i);
        }
        trid = minimal_trid();
        BCImpedance(0, 0.01, true, f, g, iPower, true, trid, file_params, modeCount);
        BCImpedance(0, 0.01, false, f, g, iPower, false, trid, file_params, modeCount);
        test.requireNear(std::real(g), 1.0, 1.0e-12, "file non-complex branch must force g");

        auto pre_params = minimal_boundary_params(BC_Mode::MODE_P_Precomputed, BC_Mode::MODE_P_Precomputed);
        pre_params.ReflectionCoef.IRC.isSet = true;
        pre_params.ReflectionCoef.IRC.xTab.resize(1);
        pre_params.ReflectionCoef.IRC.fTab.resize(1);
        pre_params.ReflectionCoef.IRC.gTab.resize(1);
        pre_params.ReflectionCoef.IRC.iTab.resize(1);
        pre_params.ReflectionCoef.IRC.xTab << 0.01;
        pre_params.ReflectionCoef.IRC.fTab << std::complex<double>(2.0, 0.0);
        pre_params.ReflectionCoef.IRC.gTab << std::complex<double>(3.0, 0.0);
        pre_params.ReflectionCoef.IRC.iTab << 7;
        trid = minimal_trid();
        BCImpedance(0, 0.01, false, f, g, iPower, false, trid, pre_params, modeCount);
        test.require(iPower == 0, "precomputed non-complex branch must reset iPower");

        Eigen::VectorXd yV(5);
        yV << 1.0, 1.0, 1.0, 1.0, 1.0;
        trid = minimal_trid();
        ElasticUP(0.01, yV, iPower, 0, trid, half_params);
        test.require(yV.allFinite(), "ElasticUP must keep finite state");
        yV << 1.0, 1.0, 1.0, 1.0, 1.0;
        ElasticDN(0.01, yV, iPower, 0, trid, half_params);
        test.require(yV.allFinite(), "ElasticDN must keep finite state");
    }

    OpenOceanKraken::ssp::SSPStructure ssp_fixture()
    {
        OpenOceanKraken::ssp::SSPStructure ssp;
        ssp.NMedia = 1;
        ssp.FirstAcoustic = 0;
        ssp.LastAcoustic = 0;
        ssp.SSPType = SSP_Mode::MODE_C_cLinear;
        ssp.NPts.resize(1);
        ssp.NMesh.resize(1);
        ssp.offset.resize(1);
        ssp.interp_offset.resize(1);
        ssp.depth.resize(1);
        ssp.beta.resize(1);
        ssp.ft.resize(1);
        ssp.sigma.resize(1);
        ssp.NPts << 4;
        ssp.NMesh << 4;
        ssp.offset << 0;
        ssp.interp_offset << 0;
        ssp.depth << 100.0;
        ssp.beta << 1.0;
        ssp.ft << 200.0;
        ssp.sigma << 0.0;
        ssp.Material = {Media_Mode::MODE_A_Acoustic};
        ssp.z.resize(4);
        ssp.alphaR.resize(4);
        ssp.alphaI.resize(4);
        ssp.betaR.resize(4);
        ssp.betaI.resize(4);
        ssp.rho.resize(4);
        ssp.z << 0.0, 25.0, 50.0, 100.0;
        ssp.alphaR << 1500.0, 1510.0, 1520.0, 1530.0;
        ssp.alphaI << 0.1, 0.2, 0.3, 0.4;
        ssp.betaR << 200.0, 220.0, 240.0, 260.0;
        ssp.betaI << 0.01, 0.02, 0.03, 0.04;
        ssp.rho << 1.0, 1.1, 1.2, 1.3;
        ssp.cp.resize(4);
        ssp.cs.resize(4);
        ssp.cpCoef = Eigen::MatrixXcd::Zero(4, 4);
        ssp.csCoef = Eigen::MatrixXcd::Zero(4, 4);
        ssp.rhoCoef = Eigen::MatrixXcd::Zero(4, 4);
        ssp.cpSpline = Eigen::MatrixXcd::Zero(4, 4);
        ssp.csSpline = Eigen::MatrixXcd::Zero(4, 4);
        ssp.rhoSpline = Eigen::MatrixXcd::Zero(4, 4);
        ssp.HSTop.BC = BC_Mode::MODE_A_Half_space;
        ssp.HSTop.alphaR = 1490.0;
        ssp.HSTop.rho = 0.9;
        ssp.HSBot.BC = BC_Mode::MODE_A_Half_space;
        ssp.HSBot.alphaR = 1600.0;
        ssp.HSBot.rho = 1.6;
        return ssp;
    }

    OOK_parameters minimal_field_params()
    {
        OOK_parameters params;
        params.freqinfo.freq = 50.0;
        params.freqinfo.Nfreq = 1;
        params.SSP = {ssp_fixture()};
        params.SourceType = Source_Mode::MODE_X_Line;
        params.coherenceType = CoherenceType::Coherent;
        params.is_Velocity = true;
        params.Pos.NSz = 1;
        params.Pos.NRz = 2;
        params.Pos.NRr = 2;
        params.Pos.NRo = 2;
        params.Pos.NRz_per_range = 2;
        params.Pos.GridType = Grid_Mode::MODE_R_Rectangular;
        params.Pos.Sz.resize(1);
        params.Pos.Rz.resize(2);
        params.Pos.Rr.resize(2);
        params.Pos.Ro.resize(2);
        params.Pos.Sz << 25.0;
        params.Pos.Rz << 0.0, 100.0;
        params.Pos.Rr << 0.0, 10.0;
        params.Pos.Ro << 0.0, 1.0;
        params.SBP.NSBPPts = 3;
        params.SBP.theta.resize(3);
        params.SBP.pat.resize(3);
        params.SBP.theta << 0.0, 45.0, 90.0;
        params.SBP.pat << 1.0, 0.5, 0.25;
        params.SBP.isSet = true;
        return params;
    }

    EigenParams minimal_eigen()
    {
        EigenParams eigen;
        eigen.M = 2;
        eigen.k.resize(2);
        eigen.PsiS.resize(2, 1);
        eigen.PsiR.resize(2, 2);
        eigen.dPsidzR.resize(2, 2);
        eigen.k << std::complex<double>(0.1, 0.01), std::complex<double>(0.2, 0.02);
        eigen.PsiS << std::complex<double>(1.0, 0.0), std::complex<double>(0.5, 0.1);
        eigen.PsiR << std::complex<double>(1.0, 0.0), std::complex<double>(0.8, 0.0),
            std::complex<double>(0.4, 0.1), std::complex<double>(0.2, 0.2);
        eigen.dPsidzR << std::complex<double>(0.1, 0.0), std::complex<double>(0.2, 0.0),
            std::complex<double>(0.3, 0.0), std::complex<double>(0.4, 0.0);
        return eigen;
    }

    void test_ssp_field_and_file_helpers(ook_test::TestRunner &test)
    {
        auto unit = [] {
            Atten_Mode mode;
            mode.attnUnit = AttenuationUnit::MODE_W_db_per_lambda;
            mode.absModel = OceanAbsorptionModel::Thorpe;
            return mode;
        }();

        auto ssp = ssp_fixture();
        TridMtx trid;
        trid.resize(8, 1, 1);
        OpenOceanKraken::ssp::UpdateSSPLoss(50.0, 50.0, ssp.NMedia, SSP_Mode::MODE_P_cPCHIP, unit, ssp);
        SSP_Mode pchip = SSP_Mode::MODE_P_cPCHIP;
        OpenOceanKraken::ssp::EvaluateSSP(trid, ssp, pchip, 0);
        test.require(std::real(trid.cp_int(0)) > 0.0, "PCHIP SSP interpolation must produce cp");

        OpenOceanKraken::ssp::UpdateSSPLoss(50.0, 50.0, ssp.NMedia, SSP_Mode::MODE_S_cCubic, unit, ssp);
        SSP_Mode cubic = SSP_Mode::MODE_S_cCubic;
        OpenOceanKraken::ssp::EvaluateSSP(trid, ssp, cubic, 0);
        test.require(std::real(trid.rho_int(0)) > 0.0, "cubic SSP interpolation must produce rho");

        SSP_Mode n2 = SSP_Mode::MODE_N_n2Linear;
        OpenOceanKraken::ssp::EvaluateSSP(trid, ssp, n2, 0);
        SSP_Mode linear = SSP_Mode::MODE_C_cLinear;
        OpenOceanKraken::ssp::EvaluateSSP(trid, ssp, linear, 0);
        SSP_Mode unknown = static_cast<SSP_Mode>(999);
        OpenOceanKraken::ssp::EvaluateSSP(trid, ssp, unknown, 0);

        auto degenerate = ssp_fixture();
        degenerate.z(1) = degenerate.z(0);
        degenerate.cp.setConstant(std::complex<double>(0.0, 1.0));
        degenerate.cs.setConstant(std::complex<double>(0.0, 1.0));
        OpenOceanKraken::ssp::EvaluateSSP(trid, degenerate, n2, 0);
        test.requireNear(std::real(trid.cp_int(0)), 1500.0, 1.0e-9, "n2Linear nonphysical N2 fallback");
        OpenOceanKraken::ssp::EvaluateSSP(trid, degenerate, linear, 0);
        test.require(std::isfinite(std::real(trid.cp_int(0))), "cLinear zero-depth guard");

        Eigen::VectorXcd cp;
        Eigen::VectorXcd cs;
        Eigen::VectorXd rho;
        OpenOceanKraken::ssp::Analytic(cp, cs, rho, 1, 4);
        test.requireNear(rho(0), 1.0, 1.0e-12, "Analytic ocean rho mismatch");
        OpenOceanKraken::ssp::Analytic(cp, cs, rho, 2, 1);
        test.requireNear(std::real(cp(0)), 1551.91, 1.0e-12, "Analytic half-space cp mismatch");
        OpenOceanKraken::ssp::Analytic(cp, cs, rho, 9, 3);
        test.requireNear(std::real(cs(0)), 2000.0, 1.0e-12, "Analytic elastic cs mismatch");

        double freq = 50.0;
        double freq0 = 50.0;
        HSInfo top;
        HSInfo bot;
        top.BC = BC_Mode::MODE_A_Half_space;
        bot.BC = BC_Mode::MODE_A_Half_space;
        top.alphaR = 1490.0;
        bot.alphaR = 1600.0;
        OpenOceanKraken::ssp::UpdateHSLoss(freq, freq0, unit, top, bot);
        test.require(std::real(top.cp) > 0.0 && std::real(bot.cp) > 0.0, "UpdateHSLoss must update half-spaces");

        OOK_parameters params = minimal_field_params();
        EigenParams eigen = minimal_eigen();
        std::vector<std::complex<float>> u(4), v(4), h(4);
        Evaluate(eigen, params, 0, 0, u.data(), v.data(), h.data());
        test.require(std::isfinite(std::real(u[0])), "coherent field output must be finite");

        params.SourceType = Source_Mode::MODE_R_Point;
        params.coherenceType = CoherenceType::Incoherent;
        std::fill(u.begin(), u.end(), std::complex<float>{});
        std::fill(v.begin(), v.end(), std::complex<float>{});
        std::fill(h.begin(), h.end(), std::complex<float>{});
        Evaluate(eigen, params, 0, 0, u.data(), v.data(), h.data());
        test.require(std::isfinite(std::real(u[1])), "incoherent field output must be finite");

        std::fill(u.begin(), u.end(), std::complex<float>{});
        field(eigen, params, u.data(), 0);
        params.SourceType = Source_Mode::MODE_X_Line;
        field(eigen, params, u.data(), 0);
        params.SourceType = static_cast<Source_Mode>(999);
        field(eigen, params, u.data(), 0);
        test.require(std::isfinite(std::real(u[3])), "legacy field helper must be finite");

        double source_rho = 0.0;
        double source_c0 = 0.0;
        auto source_ssp = ssp_fixture();
        getSourceEnv(source_ssp, -1.0, source_rho, source_c0);
        test.requireNear(source_c0, source_ssp.HSTop.alphaR, 1.0e-12, "source above top mismatch");
        getSourceEnv(source_ssp, 200.0, source_rho, source_c0);
        test.requireNear(source_c0, source_ssp.HSBot.alphaR, 1.0e-12, "source below bottom mismatch");
        getSourceEnv(source_ssp, 25.0, source_rho, source_c0);
        test.requireNear(source_rho, 1.1, 1.0e-12, "source exact layer point mismatch");
        source_ssp.NPts << 0;
        getSourceEnv(source_ssp, 30.0, source_rho, source_c0);
        test.require(source_c0 > 0.0, "source fallback nearest point mismatch");
        source_ssp.NMedia = 0;
        getSourceEnv(source_ssp, 30.0, source_rho, source_c0);
        test.requireNear(source_c0, source_ssp.HSTop.alphaR, 1.0e-12, "source empty SSP top fallback");

        const fs::path root = fs::temp_directory_path() / "ook_refcoef_file_helper";
        fs::remove_all(root);
        fs::create_directories(root);
        const fs::path env = root / "fixture.env";
        {
            std::ofstream env_file(env);
            env_file << "placeholder\n";
        }
        {
            std::ofstream brc(root / "fixture.brc");
            brc << "2\n0 0.1 0\n90 0.2 0.3\n";
        }
        OOK_parameters ref_params;
        test.require(read_refCoef_file(env.string(), ref_params, ".brc"), "read_refCoef_file .brc must succeed");
        test.require(ref_params.ReflectionCoef.RBot.size() == 2, "read_refCoef_file .brc size mismatch");
        test.require(!read_refCoef_file(env.string(), ref_params, ".unknown"), "read_refCoef_file unknown pattern must fail");
        fs::remove_all(root);
    }

    class ConcreteParamsBase : public paramsBase
    {
    public:
        void Default(OOK_parameters &) const override {}
    };

    OpenOceanKraken::ssp::SSPLayer simple_layer(int nmesh = 100)
    {
        Eigen::VectorXd z(3);
        Eigen::VectorXd alphaR(3);
        Eigen::VectorXd alphaI = Eigen::VectorXd::Zero(3);
        Eigen::VectorXd betaR = Eigen::VectorXd::Zero(3);
        Eigen::VectorXd betaI = Eigen::VectorXd::Zero(3);
        Eigen::VectorXd rho(3);
        z << 0.0, 50.0, 100.0;
        alphaR << 1500.0, 1510.0, 1520.0;
        rho << 1.0, 1.1, 1.2;
        return OpenOceanKraken::ssp::SSPLayer(3, nmesh, 0.0, 0.0, 0.0, z, alphaR, alphaI, betaR, betaI, rho, Media_Mode::MODE_A_Acoustic);
    }

    void test_input_module_helpers(ook_test::TestRunner &test)
    {
        OOK_parameters params;
        ConcreteParamsBase base;
        base.Init(params);
        base.SetupPre(params);
        base.Default(params);
        base.Read(params);
        base.Write(params);
        base.SetupPost(params);
        base.Validate(params);
        base.Echo(params);
        base.Preprocess(params);
        base.Finalize(params);

        OOK_output output{};
        outputBase out_base;
        out_base.Init(output);
        out_base.Preprocess(params, output);
        out_base.Run(params, output);
        out_base.Postprocess(params, output);
        out_base.ClearResults(params, output);
        out_base.Finalize(output);

        input_Freq freq;
        freq.Init(params);
        freq.Default(params);
        freq.Preprocess(params);
        freq.Echo(params);
        std::string title = "input module direct test";
        freq.set_title(params, title);
        Eigen::VectorXd freqvec(2);
        freqvec << 50.0, 75.0;
        freq.set_freqvec(params, freqvec);
        freq.set_freq(params, 25.0);
        freq.set_SourceType(params, Source_Mode::MODE_X_Line);
        freq.set_RunMode(params, Run_Mode::MODE_F_Field);
        freq.set_cPhase(params, 1400.0, 1700.0);
        freq.set_Velocity_enable(params, true);
        test.require(params.Title == title, "input_Freq title setter mismatch");
        test.require(params.SourceType == Source_Mode::MODE_X_Line, "input_Freq source setter mismatch");

        input_Sz_Rz_RR pos;
        pos.Init(params);
        pos.Default(params);
        Eigen::VectorXd manual(3);
        manual << 1.0, 2.0, 3.0;
        pos.set_Sz(params, manual);
        pos.set_Rz(params, manual);
        pos.set_Rr(params, manual);
        pos.set_Ro(params, manual);
        pos.set_Sz(params, 0.0, 10.0, 3);
        pos.set_Rz(params, 0.0, 20.0, 3);
        pos.set_Rr(params, 0.0, 30.0, 3);
        pos.set_Ro(params, 0.0, 40.0, 3);
        pos.set_RMax(params, 1234.0);
        pos.set_GridType(params, Grid_Mode::MODE_I_Irregular);
        pos.Preprocess(params);
        pos.Echo(params);
        test.require(params.Pos.NRz_per_range == 1, "irregular grid should set one receiver depth per range");
        test.requireNear(params.Pos.Delta_r, 15.0, 1.0e-12, "input_Sz_Rz_RR Delta_r mismatch");

        input_reflcoef refl;
        refl.Init(params);
        refl.Default(params);
        std::vector<ReflectionCoef> table = {{0.0, 0.5, 10.0}, {90.0, 0.25, 20.0}};
        refl.set_ReflCoef_Top(params, table);
        refl.set_ReflCoef_Bot(params, table);
        refl.Preprocess(params);
        auto top = refl.get_ReflCoef_Top(params);
        auto bot = refl.get_ReflCoef_Bot(params);
        refl.Preprocess(params);
        refl.Echo(params);
        test.require(top.size() == 2 && bot.size() == 2, "reflection coefficient getter size mismatch");
        test.requireNear(top[0].phi, 10.0, 1.0e-9, "reflection coefficient degree conversion mismatch");

        input_sbp sbp;
        sbp.Init(params);
        sbp.Default(params);
        Eigen::VectorXd pat(2);
        Eigen::VectorXd theta(2);
        pat << 1.0, 0.5;
        theta << -10.0, 10.0;
        sbp.set_Pat(params, pat, theta);
        sbp.Preprocess(params);
        sbp.Echo(params);
        test.require(params.SBP.isSet && params.SBP.NSBPPts == 2, "SBP setter mismatch");

        input_SSP ssp_input;
        ssp_input.Init(params);
        ssp_input.Default(params);
        ssp_input.Echo(params);
        std::vector<Range_Independent_Area> empty_areas;
        ssp_input.set_SSP(params, empty_areas);
        test.require(params.SSP.empty(), "empty SSP setter should clear SSP");

        Range_Independent_Area mutable_area;
        auto layer = simple_layer();
        mutable_area.addLayer(layer);
        mutable_area.insertLayer(1, layer);
        mutable_area.removeLayer(1);
        test.require(mutable_area.NMedia() == 1, "Range_Independent_Area layer operations mismatch");
        mutable_area.clearLayer();
        test.require(mutable_area.NMedia() == 0, "Range_Independent_Area clear mismatch");

        Range_Independent_Area area;
        area.addLayer(simple_layer());
        area.SSPType = SSP_Mode::MODE_C_cLinear;
        area.set_Top_type(BC_Mode::MODE_V_Vacuum);
        area.set_Bottom_Line(100.0, 1600.0, 0.1, 0.0, 0.0, 1.5);
        ssp_input.set_SSP(params, {area});
        Atten_Mode ssp_unit;
        ssp_unit.attnUnit = AttenuationUnit::MODE_W_db_per_lambda;
        ssp_unit.absModel = OceanAbsorptionModel::None;
        ssp_input.set_AttenUnit(params, ssp_unit);
        params.freqinfo.freq = 50.0;
        params.mesh.NSets = 2;
        params.mesh.NV[0] = 1;
        params.mesh.NV[1] = 2;
        ssp_input.Preprocess(params);
        test.require(params.NMeshMax > 0 && params.NMediaMax == 1, "input_SSP preprocess dimensions mismatch");

        auto bad_params = params;
        bad_params.SSP[0].NMesh[0] = -1;
        test.requireThrows([&]() { ssp_input.Preprocess(bad_params); }, "negative SSP mesh must throw");
    }

    void test_util_helpers(ook_test::TestRunner &test)
    {
        Eigen::VectorXd x;
        Util::linspace(2.0, 4.0, 0, x);
        Util::linspace(2.0, 4.0, 1, x);
        test.requireNear(x[0], 2.0, 1.0e-12, "linspace n=1 mismatch");
        Util::linspace(2.0, 4.0, 3, x);
        test.requireNear(x[1], 3.0, 1.0e-12, "linspace midpoint mismatch");

        Eigen::VectorXd sub(3);
        sub << 1.0, -999.9, -999.9;
        Util::SubTab(sub, 3);
        test.requireNear(sub[2], 1.0, 1.0e-12, "SubTab inherited singleton mismatch");
        Eigen::VectorXd sorted(5);
        sorted << 3.0, 1.0, 4.0, 2.0, 0.5;
        Util::Sort(sorted, 5);
        test.requireNear(sorted[0], 0.5, 1.0e-12, "Sort first element mismatch");
        Eigen::VectorXd singleton(1);
        singleton << 9.0;
        Util::Sort(singleton, 1);
        test.requireNear(singleton[0], 9.0, 1.0e-12, "Sort singleton mismatch");

        test.require(Util::spacing(1.0) > 0.0, "spacing must be positive");
        test.require(Util::isSmallValue(1.0e-5, 1.0, true), "const threshold small value");
        test.require(!Util::isSmallValue(1.0, 1.0, false), "relative threshold large value");
        auto cf = Util::cpxd2cpxf({1.25, -2.5});
        auto cd = Util::cpxf2cpxd(cf);
        test.requireNear(std::real(cd), 1.25, 1.0e-6, "complex conversion real mismatch");
        test.requireNear(std::imag(cd), -2.5, 1.0e-6, "complex conversion imag mismatch");
    }

    void test_all_env_flp_cases(ook_test::TestRunner &test, const fs::path &test_root)
    {
        std::vector<std::string> cases;
        for (const auto &entry : fs::directory_iterator(test_root))
        {
            if (entry.is_regular_file() && entry.path().extension() == ".env")
            {
                cases.push_back(entry.path().filename().string());
            }
        }
        std::sort(cases.begin(), cases.end());
        test.require(!cases.empty(), "ENV/FLP coverage requires at least one ENV case");

        int parsed = 0;
        for (const auto &name : cases)
        {
            const fs::path env_path = test_root / name;
            fs::path flp_path = env_path;
            flp_path.replace_extension(".flp");
            test.require(fs::exists(env_path), name + ": ENV file missing");
            test.require(fs::exists(flp_path), name + ": FLP file missing");

            Interface iface;
            if (is_expected_env_failure(name))
            {
                test.require(!iface.from_env((test_root / name).string()), name + ": ordinary KRAKEN should reject unsupported or missing reflection coefficient input");
                ++parsed;
                continue;
            }

            test.require(iface.from_env((test_root / name).string()), name + ": from_env must succeed");
            const auto &params = iface.getParams_const();
            require_position_consistency(test, params, name);
            require_ssp_consistency(test, params, name);
            test.require(params.runMode == Run_Mode::MODE_B_Both, name + ": FLP mode should map to Both");
            test.requireThrows([&]() { iface.get_u(0); }, name + ": get_u before run must throw");
            ++parsed;
        }
        test.require(parsed == static_cast<int>(cases.size()), "ENV/FLP coverage must include every ENV case in test directory");
    }

    void test_option_fixture_cases(ook_test::TestRunner &test, const fs::path &option_root)
    {
        if (!fs::exists(option_root))
        {
            return;
        }

        int parsed = 0;
        for (const auto &entry : fs::directory_iterator(option_root))
        {
            if (!entry.is_regular_file() || entry.path().extension() != ".env")
            {
                continue;
            }

            Interface iface;
            const auto env_path = entry.path();
            const auto name = env_path.filename().string();
            test.require(iface.from_env(env_path.string()), name + ": option fixture from_env must succeed");
            const auto &params = iface.getParams_const();
            require_position_consistency(test, params, name);
            require_ssp_consistency(test, params, name);

            const auto &area = params.sspInput.at(0);
            if (name == "option_pamf_bottom_g.env")
            {
                test.require(area.SSPType == SSP_Mode::MODE_P_cPCHIP, name + ": SSP P option mismatch");
                test.require(area.HSTop.BC == BC_Mode::MODE_A_Half_space, name + ": top A option mismatch");
                test.require(area.HSBot.BC == BC_Mode::MODE_G_Grain, name + ": bottom G option mismatch");
                test.require(params.AttenUnit.attnUnit == AttenuationUnit::MODE_M_dB_per_m, name + ": attenuation M mismatch");
                test.require(params.AttenUnit.absModel == OceanAbsorptionModel::FrancGarr, name + ": absorption F mismatch");
            }
            else if (name == "option_srnt_bottom_r.env")
            {
                test.require(area.SSPType == SSP_Mode::MODE_S_cCubic, name + ": SSP S option mismatch");
                test.require(area.HSTop.BC == BC_Mode::MODE_R_Rigid, name + ": top R option mismatch");
                test.require(area.HSBot.BC == BC_Mode::MODE_R_Rigid, name + ": bottom R option mismatch");
                test.require(params.AttenUnit.attnUnit == AttenuationUnit::MODE_N_Nepers_per_m, name + ": attenuation N mismatch");
                test.require(params.AttenUnit.absModel == OceanAbsorptionModel::Thorpe, name + ": absorption T mismatch");
            }
            else if (name == "option_cgq_bottom_v.env")
            {
                test.require(area.SSPType == SSP_Mode::MODE_C_cLinear, name + ": SSP C option mismatch");
                test.require(area.HSTop.BC == BC_Mode::MODE_G_Grain, name + ": top G option mismatch");
                test.require(area.HSBot.BC == BC_Mode::MODE_V_Vacuum, name + ": bottom V option mismatch");
                test.require(params.AttenUnit.attnUnit == AttenuationUnit::MODE_Q_Quality_Factor, name + ": attenuation Q mismatch");
            }
            else if (name == "option_nvl_bottom_p.env")
            {
                test.require(area.SSPType == SSP_Mode::MODE_N_n2Linear, name + ": SSP N option mismatch");
                test.require(area.HSBot.BC == BC_Mode::MODE_P_Precomputed, name + ": bottom P option mismatch");
                test.require(params.ReflectionCoef.IRC.isSet, name + ": IRC table must be loaded");
                test.require(params.AttenUnit.attnUnit == AttenuationUnit::MODE_L_params_lose, name + ": attenuation L mismatch");
                test.require(params.SourceType == Source_Mode::MODE_X_Line, name + ": FLP X source option mismatch");
                test.require(params.modeType == ModeType::Couple, name + ": FLP C mode option mismatch");
                test.require(params.coherenceType == CoherenceType::Incoherent, name + ": FLP I coherence option mismatch");
                test.require(!params.SBP.isSet, name + ": FLP O option should leave SBP unset");
            }
            else if (name == "option_default_chars.env")
            {
                test.require(area.SSPType == SSP_Mode::MODE_C_cLinear, name + ": unknown SSP option should default to C");
                test.require(area.HSTop.BC == BC_Mode::MODE_V_Vacuum, name + ": unknown top option should default to V");
                test.require(area.HSBot.BC == BC_Mode::MODE_V_Vacuum, name + ": unknown bottom option should keep V");
                test.require(params.AttenUnit.attnUnit == AttenuationUnit::MODE_W_db_per_lambda, name + ": unknown attenuation should default to W");
                test.require(params.AttenUnit.absModel == OceanAbsorptionModel::None, name + ": unknown absorption should default to none");
            }
            else if (name == "option_cfw_top_file.env")
            {
                test.require(area.HSTop.BC == BC_Mode::MODE_F_File, name + ": top F option mismatch");
                test.require(params.ReflectionCoef.RTop.size() == 2, name + ": top reflection table size mismatch");
            }
            else if (name == "boundary_bottom_acoustic.env")
            {
                test.require(area.HSTop.BC == BC_Mode::MODE_V_Vacuum, name + ": top V option mismatch");
                test.require(area.HSBot.BC == BC_Mode::MODE_A_Half_space, name + ": bottom A option mismatch");
                test.requireNear(area.HSBot.betaR, 0.0, 1.0e-12, name + ": bottom must be acoustic");
            }
            else if (name == "boundary_bottom_elastic.env")
            {
                test.require(area.HSBot.BC == BC_Mode::MODE_A_Half_space, name + ": bottom A option mismatch");
                test.requireNear(area.HSBot.betaR, 1100.0, 1.0e-12, name + ": bottom shear speed mismatch");
            }
            else if (name == "boundary_top_vacuum.env")
            {
                test.require(area.HSTop.BC == BC_Mode::MODE_V_Vacuum, name + ": top V option mismatch");
            }
            else if (name == "boundary_top_rigid.env")
            {
                test.require(area.HSTop.BC == BC_Mode::MODE_R_Rigid, name + ": top R option mismatch");
            }
            else if (name == "boundary_top_acoustic.env")
            {
                test.require(area.HSTop.BC == BC_Mode::MODE_A_Half_space, name + ": top A option mismatch");
                test.requireNear(area.HSTop.alphaR, 1475.0, 1.0e-12, name + ": top compressional speed mismatch");
                test.requireNear(area.HSTop.betaR, 0.0, 1.0e-12, name + ": top must be acoustic");
                test.requireNear(area.HSTop.rho, 0.98, 1.0e-12, name + ": top density mismatch");
            }
            else if (name == "boundary_top_elastic.env")
            {
                test.require(area.HSTop.BC == BC_Mode::MODE_A_Half_space, name + ": top A option mismatch");
                test.requireNear(area.HSTop.alphaR, 1475.0, 1.0e-12, name + ": top compressional speed mismatch");
                test.requireNear(area.HSTop.betaR, 100.0, 1.0e-12, name + ": top shear speed mismatch");
                test.requireNear(area.HSTop.rho, 0.98, 1.0e-12, name + ": top density mismatch");
            }
            else if (name == "boundary_top_file.env")
            {
                test.require(area.HSTop.BC == BC_Mode::MODE_F_File, name + ": top F option mismatch");
                test.require(params.ReflectionCoef.RTop.size() == 2, name + ": top reflection table size mismatch");
            }

            const fs::path json_path = fs::temp_directory_path() / ("ook_option_" + name + ".json");
            test.require(iface.to_json(json_path.string()), name + ": option fixture to_json must succeed");
            Interface json_iface;
            test.require(json_iface.from_json(json_path.string()), name + ": option fixture from_json must succeed");
            fs::remove(json_path);
            ++parsed;
        }

        test.require(parsed >= 6, "option fixture coverage requires all generated option ENV cases");
    }

    void test_kraken_reflection_option_guards(ook_test::TestRunner &test, const fs::path &test_root)
    {
        const fs::path brc_env = test_root / "neggradK_brc.env";
        if (fs::exists(brc_env))
        {
            Interface brc_iface;
            test.require(!brc_iface.from_env(brc_env.string()), "KRAKEN bottom F/.brc must fail fast");
        }

        const fs::path irc_env = test_root / "neggradK_irc.env";
        const fs::path irc_flp = test_root / "neggradK_irc.flp";
        if (!fs::exists(irc_env) || !fs::exists(irc_flp))
        {
            return;
        }

        Interface missing_irc_iface;
        test.require(!missing_irc_iface.from_env(irc_env.string()), "KRAKEN bottom P must fail when .irc is missing");

        const fs::path root = fs::temp_directory_path() / "ook_irc_fixture";
        fs::remove_all(root);
        fs::create_directories(root);
        fs::copy_file(irc_env, root / "fixture.env", fs::copy_options::overwrite_existing);
        fs::copy_file(irc_flp, root / "fixture.flp", fs::copy_options::overwrite_existing);
        {
            std::ofstream irc(root / "fixture.irc");
            irc << "'fixture irc' 500.0\n";
            irc << "3\n";
            irc << "100.0 1.0 0.5 2.0 -0.5 0\n";
            irc << "200.0 2.0 0.0 4.0 0.0 1\n";
            irc << "300.0 3.0 -0.5 6.0 0.5 2\n";
        }

        Interface irc_iface;
        test.require(irc_iface.from_env((root / "fixture.env").string()), "KRAKEN bottom P must load a valid .irc file");
        const auto &irc = irc_iface.getParams_const().ReflectionCoef.IRC;
        test.require(irc.isSet, "IRC table must be marked as set");
        test.require(irc.xTab.size() == 3, "IRC table point count mismatch");
        test.requireNear(irc.freq, 500.0, 1.0e-9, "IRC frequency mismatch");

        std::complex<double> f;
        std::complex<double> g;
        int iPower = 0;
        InterpolateIRC(std::complex<double>(50.0, 0.0), f, g, iPower, irc);
        test.requireNear(std::real(f), 1.0, 1.0e-12, "IRC lower-bound f real mismatch");
        test.requireNear(std::imag(f), 0.5, 1.0e-12, "IRC lower-bound f imag mismatch");
        test.requireNear(std::real(g), 2.0, 1.0e-12, "IRC lower-bound g real mismatch");
        test.requireNear(std::imag(g), -0.5, 1.0e-12, "IRC lower-bound g imag mismatch");
        test.require(iPower == 0, "IRC lower-bound iPower mismatch");

        InterpolateIRC(std::complex<double>(350.0, 0.0), f, g, iPower, irc);
        test.requireNear(std::real(f), 3.0, 1.0e-12, "IRC upper-bound f real mismatch");
        test.requireNear(std::imag(f), -0.5, 1.0e-12, "IRC upper-bound f imag mismatch");
        test.require(iPower == 2, "IRC upper-bound iPower mismatch");
        fs::remove_all(root);
    }

    void test_munk_exact_values(ook_test::TestRunner &test, const fs::path &test_root)
    {
        Interface iface;
        test.require(iface.from_env((test_root / "MunkK.env").string()), "MunkK from_env must succeed");
        const auto &params = iface.getParams_const();
        test.requireNear(params.freqinfo.freq, 50.0, 1.0e-9, "MunkK freq mismatch");
        test.requireNear(params.cLow, 1500.0, 1.0e-9, "MunkK cLow mismatch");
        test.requireNear(params.cHigh, 1600.0, 1.0e-9, "MunkK cHigh mismatch");
        test.requireNear(params.sspInput[0].HSBot.Depth, 5000.0, 1.0e-9, "MunkK bottom depth mismatch");
        test.requireNear(params.sspInput[0].HSBot.alphaR, 1600.0, 1.0e-9, "MunkK bottom alphaR mismatch");
        test.requireNear(params.sspInput[0].HSBot.betaR, 0.0, 1.0e-9, "MunkK bottom betaR mismatch");
        test.requireNear(params.sspInput[0].HSBot.rho, 1.8, 1.0e-9, "MunkK bottom rho mismatch");
        test.requireNear(params.sspInput[0].HSBot.alphaI, 0.8, 1.0e-9, "MunkK bottom alphaI mismatch");
        test.requireNear(params.sspInput[0].HSBot.betaI, 0.0, 1.0e-9, "MunkK bottom betaI mismatch");
        test.require(params.Pos.NSz == 1, "MunkK FLP NSz mismatch");
        test.requireNear(params.Pos.Sz[0], 1000.0, 1.0e-9, "MunkK FLP Sz mismatch");
        test.require(params.Pos.NRz == 501, "MunkK FLP NRz mismatch");
        test.requireNear(params.Pos.Rz[0], 0.0, 1.0e-9, "MunkK FLP Rz start mismatch");
        test.requireNear(params.Pos.Rz[params.Pos.NRz - 1], 5000.0, 1.0e-9, "MunkK FLP Rz end mismatch");
        test.require(params.Pos.NRr == 1001, "MunkK FLP NRr mismatch");
        test.requireNear(params.Pos.Rr[params.Pos.NRr - 1], 100000.0, 1.0e-6, "MunkK FLP Rr end mismatch");
    }

    void test_json_roundtrip_from_env(ook_test::TestRunner &test, const fs::path &test_root)
    {
        Interface env_iface;
        test.require(env_iface.from_env((test_root / "MunkK.env").string()), "MunkK from_env before JSON must succeed");
        const fs::path json_path = fs::temp_directory_path() / "ook_munk_roundtrip.json";
        test.require(env_iface.to_json(json_path.string()), "to_json must succeed for MunkK");

        Interface json_iface;
        test.require(json_iface.from_json(json_path.string()), "from_json must load MunkK snapshot");
        const auto &params = json_iface.getParams_const();
        test.requireNear(params.freqinfo.freq, 50.0, 1.0e-9, "JSON MunkK freq mismatch");
        test.require(params.SSP.size() == params.sspInput.size(), "JSON MunkK SSP conversion count mismatch");
        test.require(params.Pos.NRz == 501, "JSON MunkK NRz mismatch");
        fs::remove(json_path);
    }

    void test_json_roundtrip_from_params(ook_test::TestRunner &test)
    {
        Interface params_iface;
        std::string title = "json params interface test";
        params_iface.set_Title(title);
        params_iface.set_Freq(77.0);
        Eigen::VectorXd freqvec(2);
        freqvec << 77.0, 88.0;
        params_iface.set_freqvec(freqvec);
        params_iface.set_Sz(12.0, 24.0, 3);
        params_iface.set_Rz(0.0, 90.0, 4);
        params_iface.set_Rr(0.0, 3000.0, 4);
        params_iface.set_Ro(0.0, 30.0, 4);
        params_iface.set_GridType(Grid_Mode::MODE_I_Irregular);
        params_iface.set_SourceType(Source_Mode::MODE_X_Line);
        params_iface.set_RunMode(Run_Mode::MODE_M_Modes);
        params_iface.set_Velocity_enable(true);
        Atten_Mode atten;
        atten.attnUnit = AttenuationUnit::MODE_Q_Quality_Factor;
        atten.absModel = OceanAbsorptionModel::Thorpe;
        params_iface.set_AttenUnit(atten);
        std::vector<ReflectionCoef> refl = {{0.0, 0.5, 10.0}, {90.0, 0.25, 20.0}};
        params_iface.set_ReflCoef_Top(refl);
        params_iface.set_ReflCoef_Bottom(refl);
        Eigen::VectorXd pat(2);
        Eigen::VectorXd theta(2);
        pat << 1.0, 0.75;
        theta << -15.0, 15.0;
        params_iface.set_SBP(pat, theta);
        Range_Independent_Area area;
        area.addLayer(simple_layer());
        area.SSPType = SSP_Mode::MODE_C_cLinear;
        area.set_Top_type(BC_Mode::MODE_V_Vacuum);
        area.set_Bottom_Line(100.0, 1600.0, 0.1, 0.0, 0.0, 1.5);
        params_iface.set_SSP({area});

        const fs::path json_path = fs::temp_directory_path() / "ook_params_roundtrip.json";
        test.require(params_iface.to_json(json_path.string()), "params to_json must succeed");

        Interface json_iface;
        test.require(json_iface.from_json(json_path.string()), "params from_json must succeed");
        const auto &params = json_iface.getParams_const();
        test.require(params.Title == title, "JSON params title mismatch");
        test.require(params.freqinfo.Nfreq == 2, "JSON params frequency count mismatch");
        test.requireNear(params.freqinfo.freqvec[1], 88.0, 1.0e-9, "JSON params freqvec mismatch");
        test.require(params.Pos.NSz == 3, "JSON params NSz mismatch");
        test.require(params.Pos.NRz == 4, "JSON params NRz mismatch");
        test.require(params.Pos.NRr == 4, "JSON params NRr mismatch");
        test.require(params.Pos.GridType == Grid_Mode::MODE_I_Irregular, "JSON params grid mode mismatch");
        test.require(params.SourceType == Source_Mode::MODE_X_Line, "JSON params source type mismatch");
        test.require(params.runMode == Run_Mode::MODE_M_Modes, "JSON params run mode mismatch");
        test.require(params.is_Velocity, "JSON params velocity flag mismatch");
        test.require(params.SBP.isSet, "JSON params SBP mismatch");
        test.require(params.SSP.size() == params.sspInput.size(), "JSON params SSP conversion count mismatch");

        const auto base_json = OpenOcean_json::parse(params_iface.to_json_string());
        auto expect_bad_json = [&](OpenOcean_json bad, const std::string &suffix) {
            const fs::path bad_path = fs::temp_directory_path() / ("ook_bad_json_" + suffix + ".json");
            {
                std::ofstream bad_file(bad_path);
                bad_file << bad.dump(2);
            }
            Interface bad_iface;
            test.require(!bad_iface.from_json(bad_path.string()), "bad JSON variant must fail: " + suffix);
            fs::remove(bad_path);
        };

        auto bad = base_json;
        bad["AttenUnit"]["AttenuationUnit"] = 12;
        expect_bad_json(bad, "atten_unit_type");
        bad = base_json;
        bad["AttenUnit"]["OceanAbsorptionModel"] = 12;
        expect_bad_json(bad, "abs_model_type");
        bad = base_json;
        bad["SourceType"] = 12;
        expect_bad_json(bad, "source_type");
        bad = base_json;
        bad["RunMode"] = 12;
        expect_bad_json(bad, "run_mode_type");
        bad = base_json;
        bad["CoherenceType"] = 12;
        expect_bad_json(bad, "coherence_type");
        bad = base_json;
        bad["ModeType"] = 12;
        expect_bad_json(bad, "mode_type");
        bad = base_json;
        bad["sspInput"][0]["layers"][0]["Material"] = 12;
        expect_bad_json(bad, "media_type");
        fs::remove(json_path);
    }

    void test_bad_input_cases(ook_test::TestRunner &test, const fs::path &bad_root)
    {
        if (!fs::exists(bad_root))
        {
            return;
        }

        int env_cases = 0;
        int json_cases = 0;
        for (const auto &entry : fs::directory_iterator(bad_root))
        {
            if (!entry.is_regular_file())
            {
                continue;
            }

            const fs::path path = entry.path();
            const std::string name = path.filename().string();
            if (path.extension() == ".env")
            {
                Interface iface;
                test.require(!iface.from_env(path.string()), name + ": bad ENV fixture must fail from_env");
                ++env_cases;
            }
            else if (path.extension() == ".json" && name != "_manifest.json")
            {
                Interface iface;
                test.require(!iface.from_json(path.string()), name + ": bad JSON fixture must fail from_json");
                ++json_cases;
            }
        }

        test.require(env_cases > 0, "LCOV bad input coverage requires at least one bad ENV fixture");
        test.require(json_cases > 0, "LCOV bad input coverage requires at least one bad JSON fixture");
    }

    void test_multilayer_elastic_stack_terminal_mode(ook_test::TestRunner &test, const fs::path &test_root)
    {
        ThreadPool pool(1);
        Interface iface(pool);
        iface.setNumThreads(1);

        const fs::path env_path = test_root / "multilayer_elastic_stack.env";
        test.require(fs::exists(env_path), "multilayer elastic stack fixture must exist");
        test.require(iface.from_env(env_path.string()), "multilayer elastic stack from_env must succeed");
        iface.runEigen();

        const auto &eigen = iface.getOutput_const().eigen[0];
        test.require(eigen.M == 19, "multilayer elastic stack must retain the terminal cutoff mode");

        const double omega = 2.0 * std::acos(-1.0) * iface.getParams_const().freqinfo.freq;
        const double terminal_phase_speed = omega / std::real(eigen.k(eigen.M - 1));
        test.requireNear(terminal_phase_speed, 1998.731201, 0.05,
                         "multilayer elastic stack terminal mode phase speed mismatch");
    }

    void test_multilayer_mud_sand_lossy_wavenumbers(ook_test::TestRunner &test, const fs::path &test_root)
    {
        ThreadPool pool(1);
        Interface iface(pool);
        iface.setNumThreads(1);

        const fs::path env_path = test_root / "multilayer_mud_sand.env";
        test.require(fs::exists(env_path), "multilayer mud/sand fixture must exist");
        test.require(iface.from_env(env_path.string()), "multilayer mud/sand from_env must succeed");
        iface.runEigen();

        const auto &eigen = iface.getOutput_const().eigen[0];
        test.require(eigen.M == 7, "multilayer mud/sand must retain seven propagating modes");

        const double kraken_k[] = {
            0.418148130178,
            0.415981858969,
            0.412391990423,
            0.406430929899,
            0.405105829239,
            0.401739627123,
            0.396181911230,
        };
        for (int mode = 0; mode < eigen.M; ++mode)
        {
            test.requireNear(std::real(eigen.k(mode)), kraken_k[mode], 1.0e-7,
                             "multilayer mud/sand lossy wavenumber mismatch at mode " + std::to_string(mode + 1));
        }
    }

    void test_mod_export_header(ook_test::TestRunner &test, const fs::path &test_root)
    {
        ThreadPool pool(1);
        Interface iface(pool);
        iface.setNumThreads(1);
        test.require(iface.from_env((test_root / "MunkK.env").string()), "MunkK from_env before MOD export must succeed");
        iface.runEigen();

        const fs::path root = fs::temp_directory_path() / "ook_munk_mod_export";
        const fs::path mod_path = root.string() + ".mod";
        fs::remove(mod_path);
        iface.export_mod(root.string());
        test.require(fs::exists(mod_path), "export_mod must create a .mod file");

        std::ifstream mod(mod_path, std::ios::binary);
        test.require(mod.is_open(), "exported .mod must be readable");

        int lrecl = 0;
        mod.read(reinterpret_cast<char *>(&lrecl), sizeof(lrecl));
        test.require(lrecl >= 32, "MOD LRecordLength must respect Kraken minimum");
        const std::streamoff rec_bytes = static_cast<std::streamoff>(4 * lrecl);

        std::vector<char> header(static_cast<size_t>(rec_bytes), 0);
        mod.seekg(0, std::ios::beg);
        mod.read(header.data(), static_cast<std::streamsize>(header.size()));
        int nfreq = 0, nmedia = 0, ntot = 0, nmat = 0;
        std::memcpy(&nfreq, header.data() + 84, sizeof(nfreq));
        std::memcpy(&nmedia, header.data() + 88, sizeof(nmedia));
        std::memcpy(&ntot, header.data() + 92, sizeof(ntot));
        std::memcpy(&nmat, header.data() + 96, sizeof(nmat));
        test.require(nfreq == 1, "OOK MOD export should write one solved frequency");
        test.require(nmedia == 1, "MunkK MOD NMedia mismatch");
        test.require(ntot == 1001, "MunkK MOD zTab size mismatch");
        test.require(nmat == ntot, "KRAKEN MOD NTot/NMat should match for acoustic modes");

        int modes = 0;
        mod.seekg(5 * rec_bytes, std::ios::beg);
        mod.read(reinterpret_cast<char *>(&modes), sizeof(modes));
        test.require(modes > 0, "MunkK MOD must contain at least one mode");
        test.require(fs::file_size(mod_path) >= static_cast<uintmax_t>((8 + modes) * rec_bytes), "MOD file is too small for mode and wavenumber records");
        mod.close();
        fs::remove(mod_path);
    }
}

int main(int argc, char **argv)
{
    ook_test::TestRunner test;
    try
    {
        const fs::path test_root = (argc > 1) ? fs::path(argv[1]) : fs::path("test");
        test.require(fs::exists(test_root), "test root must exist");
        const fs::path env_root = positive_fixture_root(test_root);
        const fs::path bad_root = test_root / "bad_examples";
        const fs::path option_root = test_root / "option_examples";
        test.require(fs::exists(env_root), "positive ENV fixture root must exist");
        test_lifecycle_and_guards(test);
        test_top_bottom_line(test);
        test_interpolation_helpers(test);
        test_attenuation_and_reflection_helpers(test);
        test_boundary_and_scatter_helpers(test);
        test_ssp_field_and_file_helpers(test);
        test_input_module_helpers(test);
        test_util_helpers(test);
        test_all_env_flp_cases(test, env_root);
        test_option_fixture_cases(test, option_root);
        test_kraken_reflection_option_guards(test, env_root);
        test_munk_exact_values(test, env_root);
        test_json_roundtrip_from_env(test, env_root);
        test_json_roundtrip_from_params(test);
        test_bad_input_cases(test, bad_root);
        test_multilayer_elastic_stack_terminal_mode(test, env_root);
        test_multilayer_mud_sand_lossy_wavenumbers(test, env_root);
        test_mod_export_header(test, env_root);
        std::cout << "OOK interface tests passed. checks=" << test.checks() << std::endl;
        std::cout << "ENV/FLP coverage=all discovered ENV cases, JSON coverage=2/2, bad fixtures=all discovered ENV/JSON" << std::endl;
        return 0;
    }
    catch (const std::exception &e)
    {
        std::cerr << "OOK interface tests failed: " << e.what() << std::endl;
        std::cerr << "checks before failure=" << test.checks() << std::endl;
        return 1;
    }
}
