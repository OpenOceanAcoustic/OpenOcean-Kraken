#include "test.h"
#include "AttenMod.h"
#include "BCImpedanceMod.h"
#include "OpenOceanKrakenKernelInterface.h"
#include "pchipMod.h"
#include "RefCoef.h"
#include "run.h"
#include "Scatter.h"
#include "splinec.h"
#include "sspMod.h"
#include "ThreadPool.h"
#include "util.h"
#include "field.h"
#include "EvaluateAD.h"
#include "EvaluateCM.h"
#include "MergeVectors.h"
#include "env_in_out.hpp"
#include "input_Freq.hpp"
#include "input_reflcoef.hpp"
#include "input_sbp.hpp"
#include "input_SSP.hpp"
#include "input_Sz_Rz_RR.hpp"
#include "output_eigen.hpp"
#include "output_field.hpp"

#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <vector>

namespace OpenOceanKraken
{
    bool read_refCoef_file(const std::string &envPath, OOK_parameters &params, std::string pattern);
}

namespace fs = std::filesystem;
using OpenOceanKraken::Atten_Mode;
using OpenOceanKraken::AttenuationUnit;
using OpenOceanKraken::BiologicalAttenuationLayer;
using OpenOceanKraken::BCImpedance;
using OpenOceanKraken::BC_Mode;
using OpenOceanKraken::CoherenceType;
using OpenOceanKraken::CRCI;
using OpenOceanKraken::CSpline;
using OpenOceanKraken::ElasticDN;
using OpenOceanKraken::ElasticUP;
using OpenOceanKraken::EigenParams;
using OpenOceanKraken::Evaluate;
using OpenOceanKraken::EvaluateAD;
using OpenOceanKraken::EvaluateCM;
using OpenOceanKraken::Grid_Mode;
using OpenOceanKraken::HSInfo;
using OpenOceanKraken::KernelInterface;
using OpenOceanKraken::InterpolateIRC;
using OpenOceanKraken::InterpolateReflectionCoefficient;
using OpenOceanKraken::InternalReflectionCoefInfo;
using OpenOceanKraken::KupIng;
using OpenOceanKraken::Media_Mode;
using OpenOceanKraken::ModeType;
using OpenOceanKraken::MaxBioLayers;
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
using OpenOceanKraken::attenuationModesEqual;
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
using OpenOceanKraken::read_env_file;
using OpenOceanKraken::read_refCoef_file;
using OpenOceanKraken::spline;
using OpenOceanKraken::ssp::Range_Independent_Area;
using OpenOceanKraken::validateAttenuationMode;
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

    std::string biological_env_text(
        const std::string &top_option,
        const std::string &biological_block)
    {
        std::ostringstream out;
        out << "'Biological ENV test'\n"
            << "1000.0\n"
            << "1\n"
            << "'" << top_option << "'\n"
            << biological_block;
        if (top_option.size() > 1 && top_option[1] == 'A')
            out << "0.0 1450.0 0.0 1.0 0.0 0.0 /\n";
        out << "1000 0.0 100.0\n"
            << "0.0 1500.0 /\n"
            << "100.0 1500.0 /\n"
            << "'A' 0.0\n"
            << "100.0 1600.0 0.0 1.5 0.0 /\n"
            << "1400.0 20000.0\n"
            << "10.0\n"
            << "1\n"
            << "50.0 /\n"
            << "1\n"
            << "50.0 /\n";
        return out.str();
    }

    const char *single_profile_flp =
        "/,\n"
        "'RA'\n"
        "1,\n"
        "1\n"
        "0.0 /\n"
        "3\n"
        "1.0 5.0 10.0 /\n"
        "1\n"
        "50.0 /\n"
        "1\n"
        "50.0 /\n"
        "1\n"
        "0.0 /\n";

    void write_text(const fs::path &path, const std::string &text)
    {
        std::ofstream stream(path);
        stream << text;
        if (!stream)
            throw std::runtime_error("cannot write test fixture " + path.string());
    }

    class ScopedTemporaryDirectory
    {
    public:
        explicit ScopedTemporaryDirectory(fs::path path)
            : path_(std::move(path))
        {
            fs::remove_all(path_);
            fs::create_directories(path_);
        }

        ~ScopedTemporaryDirectory()
        {
            std::error_code ignored;
            fs::remove_all(path_, ignored);
        }

        const fs::path &path() const
        {
            return path_;
        }

    private:
        fs::path path_;
    };

    void test_biological_env_transport(ook_test::TestRunner &test)
    {
        ScopedTemporaryDirectory root(
            fs::temp_directory_path() / "ook_biological_env_transport");

        struct EnvVariant
        {
            std::string name;
            std::string option;
            std::string block;
            bool expected;
        };

        const std::vector<EnvVariant> variants = {
            {"normal_top", "CVWB", "1\n20 40 1000 5 0.04\n", true},
            {"top_halfspace", "CAWB", "1\n20 40 1000 5 0.04\n", true},
            {"empty_layers", "CVWB", "0\n", true},
            {"missing_count", "CVWB", "", false},
            {"negative_count", "CVWB", "-1\n", false},
            {"too_many", "CVWB", "201\n", false},
            {"count_mismatch", "CVWB", "2\n20 40 1000 5 0.04\n", false},
            {"four_fields", "CVWB", "1\n20 40 1000 5\n", false},
            {"six_fields", "CVWB", "1\n20 40 1000 5 0.04 9\n", false},
            {"non_numeric", "CVWB", "1\n20 40 bad 5 0.04\n", false},
            {"invalid_depths", "CVWB", "1\n40 20 1000 5 0.04\n", false},
            {"unknown_option", "CVWX", "", false},
        };

        for (const auto &variant : variants)
        {
            const fs::path env_path = root.path() / (variant.name + ".env");
            fs::path flp_path = env_path;
            flp_path.replace_extension(".flp");
            write_text(env_path, biological_env_text(variant.option, variant.block));
            write_text(flp_path, single_profile_flp);

            KernelInterface iface;
            if (variant.expected)
            {
                test.require(iface.from_env(env_path.string()),
                             variant.name + ": Biological ENV must load");
                const auto &mode = iface.getParams_const().AttenUnit;
                test.require(mode.attnUnit == AttenuationUnit::MODE_W_db_per_lambda,
                             variant.name + ": attenuation unit mismatch");
                test.require(mode.absModel == OceanAbsorptionModel::Biological,
                             variant.name + ": Biological model mismatch");
                if (variant.block == "0\n")
                {
                    test.require(mode.biologicalLayers.empty(),
                                 variant.name + ": empty Biological layer set mismatch");
                }
                else
                {
                    test.require(mode.biologicalLayers.size() == 1,
                                 variant.name + ": Biological layer count mismatch");
                    const auto &layer = mode.biologicalLayers.front();
                    test.requireNear(layer.Z1, 20.0, 1.0e-12,
                                     variant.name + ": Biological Z1 mismatch");
                    test.requireNear(layer.Z2, 40.0, 1.0e-12,
                                     variant.name + ": Biological Z2 mismatch");
                    test.requireNear(layer.f0, 1000.0, 1.0e-12,
                                     variant.name + ": Biological f0 mismatch");
                    test.requireNear(layer.Q, 5.0, 1.0e-12,
                                     variant.name + ": Biological Q mismatch");
                    test.requireNear(layer.a0, 0.04, 1.0e-12,
                                     variant.name + ": Biological a0 mismatch");
                }
                if (variant.option[1] == 'A')
                {
                    test.requireNear(iface.getParams_const().sspInput[0].HSTop.alphaR,
                                     1450.0, 1.0e-12,
                                     variant.name + ": top half-space must follow Biological block");
                }
            }
            else
            {
                Atten_Mode thorpe;
                thorpe.absModel = OceanAbsorptionModel::Thorpe;
                iface.set_AttenUnit(thorpe);
                test.require(!iface.from_env(env_path.string()),
                             variant.name + ": invalid Biological ENV must fail");
                test.require(attenuationModesEqual(
                                 iface.getParams_const().AttenUnit, thorpe),
                             variant.name + ": failed ENV import must preserve attenuation mode");
            }
        }

        const std::string one_layer = "1\n20 40 1000 5 0.04\n";
        const std::string two_layers =
            "2\n20 40 1000 5 0.04\n50 80 2000 6 0.05\n";
        const fs::path profiles_path = root.path() / "profiles.env";
        write_text(profiles_path,
                   biological_env_text("CVWB", one_layer) +
                       biological_env_text("CVWB", one_layer));
        OOK_parameters matching_profiles;
        test.require(read_env_file(profiles_path.string(), matching_profiles),
                     "identical Biological ENV profiles must load");
        test.require(matching_profiles.sspInput.size() == 2,
                     "identical Biological ENV profiles must retain both SSPs");
        test.require(matching_profiles.AttenUnit.absModel ==
                         OceanAbsorptionModel::Biological,
                     "identical Biological ENV profiles must retain Biological mode");

        struct InconsistentProfiles
        {
            std::string name;
            std::string first_option;
            std::string first_block;
            std::string second_option;
            std::string second_block;
        };

        const std::vector<InconsistentProfiles> inconsistent_profiles = {
            {"unit", "CVWB", one_layer, "CVMB", one_layer},
            {"model", "CVWB", one_layer, "CVWT", ""},
            {"layer_count", "CVWB", one_layer, "CVWB", two_layers},
            {"layer_field", "CVWB", two_layers, "CVWB",
             "2\n20 40 1000 5 0.04\n50 80 2000 6 0.06\n"},
            {"layer_order", "CVWB", two_layers, "CVWB",
             "2\n50 80 2000 6 0.05\n20 40 1000 5 0.04\n"},
        };

        for (const auto &variant : inconsistent_profiles)
        {
            const fs::path env_path = root.path() /
                                      ("inconsistent_" + variant.name + ".env");
            write_text(env_path,
                       biological_env_text(variant.first_option, variant.first_block) +
                           biological_env_text(variant.second_option, variant.second_block));
            OOK_parameters parsed;
            test.require(!read_env_file(env_path.string(), parsed),
                         variant.name + ": inconsistent attenuation profiles must fail");
        }
    }

    void test_singleton_mode_position_interpolation(
        ook_test::TestRunner &test)
    {
        ScopedTemporaryDirectory root(
            fs::temp_directory_path() / "ook_singleton_mode_position");
        const fs::path env_path = root.path() / "singleton.env";
        fs::path flp_path = env_path;
        flp_path.replace_extension(".flp");
        write_text(env_path, biological_env_text("CVW", ""));
        write_text(flp_path, single_profile_flp);

        ThreadPool pool(1);
        KernelInterface iface(pool);
        iface.setNumThreads(1);
        test.require(
            iface.from_env(env_path.string()),
            "singleton ModePos ENV must load");

        const auto &params = iface.getParams_const();
        test.require(params.hasModePos, "ENV load must preserve ModePos");
        test.require(
            params.ModePos.Sz.size() == 1 &&
                params.ModePos.Rz.size() == 1,
            "singleton ModePos must contain one source and receiver depth");
        test.requireNear(
            params.ModePos.Sz(0), 50.0, 0.0,
            "singleton ModePos source depth mismatch");
        test.requireNear(
            params.ModePos.Rz(0), 50.0, 0.0,
            "singleton ModePos receiver depth mismatch");

        iface.runEigen();
        const auto &eigen = iface.getOutput_const().eigen[0];
        test.require(eigen.M > 0, "singleton ModePos must produce modes");
        test.require(
            eigen.ModeZ.size() == 1 && eigen.PhiMode.cols() == 1,
            "singleton ModePos must retain one interpolation column");
        for (int mode = 0; mode < eigen.M; ++mode)
        {
            const auto psi_s = eigen.PsiS(mode, 0);
            const auto psi_r = eigen.PsiR(mode, 0);
            const auto phi = eigen.PhiMode(mode, 0);
            const auto dpsi_s = eigen.dPsidzS(mode, 0);
            const auto dpsi_r = eigen.dPsidzR(mode, 0);
            test.require(
                std::isfinite(psi_s.real()) && std::isfinite(psi_s.imag()) &&
                    std::isfinite(psi_r.real()) &&
                    std::isfinite(psi_r.imag()),
                "singleton ModePos eigenfunctions must be finite");
            test.requireNear(
                std::abs(psi_s - phi), 0.0, 1.0e-12,
                "singleton source interpolation must use PhiMode col0");
            test.requireNear(
                std::abs(psi_r - phi), 0.0, 1.0e-12,
                "singleton receiver interpolation must use PhiMode col0");
            test.require(
                std::isfinite(dpsi_s.real()) &&
                    std::isfinite(dpsi_s.imag()) &&
                    std::isfinite(dpsi_r.real()) &&
                    std::isfinite(dpsi_r.imag()),
                "singleton ModePos derivatives must be finite");
            test.requireNear(
                std::abs(dpsi_s - dpsi_r), 0.0, 1.0e-12,
                "singleton source and receiver derivatives must agree");
        }

        iface.runField();
        const auto pressure = iface.getPressureCopy();
        test.require(
            pressure.values.size() == 3,
            "singleton fixture must retain three field ranges");
        for (const auto value : pressure.values)
        {
            test.require(
                std::isfinite(value.real()) && std::isfinite(value.imag()),
                "singleton ModePos field must be finite");
        }
    }

    void test_lifecycle_and_guards(ook_test::TestRunner &test)
    {
        KernelInterface iface;
        test.require(iface.getNumThreads() == 1, "Default KernelInterface must use one thread");
        test.requireThrows([&]() { iface.setNumThreads(0); }, "setNumThreads(0) must throw");
        test.requireThrows([&]() { iface.runEigen(); }, "runEigen without ThreadPool must throw");
        test.requireThrows([&]() { iface.getOutput_Copy(); }, "getOutput_Copy must be disabled");

        KernelInterface freed;
        freed.free();
        freed.free();
        test.requireThrows([&]() { freed.getParams(); }, "getParams after free must throw");
    }

    void test_interface_velocity_and_export_workflow(ook_test::TestRunner &test, const fs::path &test_root)
    {
        const fs::path env_path = test_root / "calibK.env";
        if (!fs::exists(env_path))
        {
            return;
        }

        ThreadPool pool(1);
        KernelInterface iface(pool);
        iface.setThreadPool(pool);
        iface.setNumThreads(1);
        test.require(iface.getHardwareThreads() >= 0, "hardware thread query must be available");
        test.require(iface.from_env(env_path.string()), "interface velocity fixture must load");

        std::string title = "interface velocity and export workflow";
        iface.set_Title(title);
        iface.set_Freq(250.0);
        Eigen::VectorXd freqvec(1);
        freqvec << 250.0;
        iface.set_freqvec(freqvec);
        Eigen::VectorXd sz(1);
        Eigen::VectorXd rz(3);
        Eigen::VectorXd rr(3);
        Eigen::VectorXd ro(3);
        sz << 50.0;
        rz << 0.0, 50.0, 100.0;
        rr << 1.0, 2.0, 3.0;
        ro.setZero();
        iface.set_Sz(sz);
        iface.set_Rz(rz);
        iface.set_Rr(rr);
        iface.set_Ro(ro);
        iface.set_cPhase(1400.0, 20000.0);
        iface.set_Rmax(1000.0);
        iface.set_GridType(Grid_Mode::MODE_R_Rectangular);
        iface.set_SourceType(Source_Mode::MODE_X_Line);
        iface.set_RunMode(Run_Mode::MODE_B_Both);
        iface.set_Velocity_enable(true);
        iface.run();

        test.require(iface.get_u(0) != nullptr, "source pressure getter must return data");
        test.require(iface.get_v(0) != nullptr, "source vertical velocity getter must return data");
        test.require(iface.get_h(0) != nullptr, "source horizontal velocity getter must return data");
        test.require(iface.get_u_AllSources() != nullptr, "all-source pressure getter must return data");
        test.require(iface.get_v_AllSources() != nullptr, "all-source vertical velocity getter must return data");
        test.require(iface.get_h_AllSources() != nullptr, "all-source horizontal velocity getter must return data");
        test.requireThrows([&]() { iface.get_u(1); }, "source getter must reject an out-of-range source");

        const fs::path root = fs::temp_directory_path() / "ook_interface_velocity_export";
        fs::remove_all(root);
        fs::create_directories(root);
        const std::string output_root = (root / "result").string();
        iface.export_result(output_root);
        iface.export_shd((root / "default_component").string(), 0);
        test.require(fs::exists(root / "result_P.mod"), "export_result must write pressure MOD");
        test.require(fs::exists(root / "result_P.shd"), "export_result must write pressure SHD");
        test.require(fs::exists(root / "result_V.shd"), "export_result must write vertical velocity SHD");
        test.require(fs::exists(root / "result_H.shd"), "export_result must write horizontal velocity SHD");
        test.require(fs::exists(root / "default_component.shd"), "default SHD component must write pressure");

        iface.clearResults();
        test.requireNear(std::abs(iface.get_u(0)[0]), 0.0, 0.0, "clearResults must zero pressure");
        test.requireNear(std::abs(iface.get_v(0)[0]), 0.0, 0.0, "clearResults must zero vertical velocity");
        test.requireNear(std::abs(iface.get_h(0)[0]), 0.0, 0.0, "clearResults must zero horizontal velocity");
        fs::remove_all(root);
    }

    void test_top_bottom_line(ook_test::TestRunner &test)
    {
        Range_Independent_Area area;
        area.set_Bottom_Line(5000.0, 1600.0, 0.8, 0.0, 0.0, 1.8);
        area.set_Top_Line(0.0, 1500.0, 0.0, 0.0, 0.0, 1.0);
        area.set_Bottom_type(BC_Mode::MODE_R_Rigid);
        area.set_Bottom_type(BC_Mode::MODE_A_Half_space);
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
        test.require(addOceanAbsorption(0.0, 0.0, 1000.0, unit(AttenuationUnit::MODE_W_db_per_lambda, OceanAbsorptionModel::Thorpe)) > 0.0, "Thorpe absorption must be positive");
        test.require(addOceanAbsorption(0.0, 0.0, 1000.0, unit(AttenuationUnit::MODE_W_db_per_lambda, OceanAbsorptionModel::FrancGarr)) > 0.0, "Franc-Garr absorption must be positive");

        Atten_Mode biological;
        biological.absModel = OceanAbsorptionModel::Biological;
        biological.biologicalLayers = {{20.0, 40.0, 1000.0, 5.0, 0.04}};

        constexpr double expected_np_per_m = 1.0 / 8685.8896;
        test.requireNear(addOceanAbsorption(0.0, 20.0, 1000.0, biological), expected_np_per_m, 1.0e-15, "lower boundary must be included");
        test.requireNear(addOceanAbsorption(0.0, 30.0, 1000.0, biological), expected_np_per_m, 1.0e-15, "layer interior mismatch");
        test.requireNear(addOceanAbsorption(0.0, 40.0, 1000.0, biological), expected_np_per_m, 1.0e-15, "upper boundary must be included");
        test.requireNear(addOceanAbsorption(0.0, 30.0, 500.0, biological), 5.0942148298338349e-7, 1.0e-18, "below-resonance f=500 Biological attenuation must use the squared frequency ratio");
        test.requireNear(addOceanAbsorption(0.0, 30.0, 2000.0, biological), 7.64343602683782e-6, 1.0e-18, "above-resonance f=2000 Biological attenuation must use the squared frequency ratio");
        test.requireNear(addOceanAbsorption(0.0, 19.999, 1000.0, biological), 0.0, 1.0e-15, "point below layer must be excluded");
        test.requireNear(addOceanAbsorption(0.0, 40.001, 1000.0, biological), 0.0, 1.0e-15, "point above layer must be excluded");

        Atten_Mode empty_biological = biological;
        empty_biological.biologicalLayers.clear();
        test.requireNear(addOceanAbsorption(0.25, 30.0, 1000.0, empty_biological), 0.25, 0.0, "empty Biological layers must add no attenuation");

        Atten_Mode zero_thickness = biological;
        zero_thickness.biologicalLayers = {{30.0, 30.0, 1000.0, 5.0, 0.04}, {20.0, 40.0, 1000.0, 5.0, 0.0}};
        test.requireNear(addOceanAbsorption(0.0, 30.0, 1000.0, zero_thickness), expected_np_per_m, 1.0e-15, "zero-thickness layers are valid and a0 zero adds nothing");

        Atten_Mode overlapping = biological;
        overlapping.biologicalLayers.push_back({30.0, 50.0, 1000.0, 5.0, 0.04});
        test.requireNear(addOceanAbsorption(0.25, 30.0, 1000.0, overlapping), 0.25 + 2.0 * expected_np_per_m, 1.0e-15, "overlapping layers and shared endpoints must accumulate");

        double z = 0.0;
        double c = 1500.0;
        double alpha = 0.1;
        double freq = 100.0;
        double freq0 = 100.0;
        double beta = 1.0;
        double ft = 200.0;
        auto crci = CRCI(z, c, alpha, freq, freq0, unit(AttenuationUnit::MODE_W_db_per_lambda), beta, ft);
        test.requireNear(std::real(crci), 1500.0, 1.0e-12, "CRCI real component mismatch");

        double bio_z = 30.0;
        double bio_c = 1500.0;
        double bio_alpha = 0.0;
        double bio_freq = 1000.0;
        double bio_freq0 = 1000.0;
        double bio_beta = 1.0;
        double bio_ft = 2000.0;
        const auto bio_crci = CRCI(bio_z, bio_c, bio_alpha, bio_freq, bio_freq0, biological, bio_beta, bio_ft);
        test.requireNear(std::imag(bio_crci), 0.041227627617643738, 1.0e-12, "weak-loss Biological complex sound speed mismatch");
        test.requireThrows([&] { double zero_frequency = 0.0; CRCI(bio_z, bio_c, bio_alpha, zero_frequency, bio_freq0, biological, bio_beta, bio_ft); }, "zero attenuation frequency must fail");

        Atten_Mode overflowing = biological;
        overflowing.biologicalLayers[0].f0 = std::numeric_limits<double>::max();
        test.requireThrows([&] { addOceanAbsorption(0.0, 30.0, 1000.0, overflowing); }, "non-finite Biological denominator must fail");
        overflowing = biological;
        overflowing.biologicalLayers[0].a0 = std::numeric_limits<double>::max();
        test.requireThrows([&] { addOceanAbsorption(0.0, 30.0, 1000.0, overflowing); }, "non-finite Biological contribution must fail");

        Atten_Mode fatal_loss = biological;
        fatal_loss.biologicalLayers[0].a0 = 1.0e8;
        test.requireThrows([&] { CRCI(bio_z, bio_c, bio_alpha, bio_freq, bio_freq0, fatal_loss, bio_beta, bio_ft); }, "imaginary sound speed greater than real sound speed must fail");
        double fluid_shear_speed = 0.0;
        const auto fluid_shear = CRCI(bio_z, fluid_shear_speed, bio_alpha, bio_freq, bio_freq0, biological, bio_beta, bio_ft);
        test.requireNear(std::abs(fluid_shear), 0.0, 0.0, "fluid cS zero must remain valid");

        KernelInterface biological_instance;
        KernelInterface default_instance;
        biological_instance.set_AttenUnit(biological);
        test.require(biological_instance.getParams_const().AttenUnit.absModel == OceanAbsorptionModel::Biological, "first instance must retain Biological mode");
        test.require(default_instance.getParams_const().AttenUnit.absModel == OceanAbsorptionModel::None && default_instance.getParams_const().AttenUnit.biologicalLayers.empty(), "Biological state must not leak between KernelInterface instances");

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

    void test_biological_attenuation_validation(ook_test::TestRunner &test)
    {
        Atten_Mode valid;
        valid.absModel = OceanAbsorptionModel::Biological;
        valid.biologicalLayers = {{20.0, 40.0, 1000.0, 5.0, 0.04}};
        validateAttenuationMode(valid);

        Atten_Mode empty = valid;
        empty.biologicalLayers.clear();
        validateAttenuationMode(empty);

        Atten_Mode copy = valid;
        test.require(attenuationModesEqual(valid, copy),
                     "identical Biological modes must compare equal");
        std::swap(copy.biologicalLayers[0].Z1, copy.biologicalLayers[0].Z2);
        test.require(!attenuationModesEqual(valid, copy),
                     "layer field or order changes must compare unequal");

        auto require_invalid = [&](Atten_Mode candidate, const std::string &name) {
            test.requireThrows(
                [&] { validateAttenuationMode(candidate); },
                name + " must be rejected");
        };

        Atten_Mode candidate = valid;
        candidate.biologicalLayers[0].Z1 =
            std::numeric_limits<double>::quiet_NaN();
        require_invalid(candidate, "non-finite Z1");
        candidate = valid;
        candidate.biologicalLayers[0].Z1 = 41.0;
        require_invalid(candidate, "Z1 greater than Z2");
        candidate = valid;
        candidate.biologicalLayers[0].f0 = 0.0;
        require_invalid(candidate, "non-positive f0");
        candidate = valid;
        candidate.biologicalLayers[0].Q = 0.0;
        require_invalid(candidate, "non-positive Q");
        candidate = valid;
        candidate.biologicalLayers[0].a0 = -0.01;
        require_invalid(candidate, "negative a0");
        candidate = valid;
        candidate.biologicalLayers.resize(MaxBioLayers + 1, valid.biologicalLayers[0]);
        require_invalid(candidate, "more than 200 layers");
        candidate = valid;
        candidate.absModel = OceanAbsorptionModel::Thorpe;
        require_invalid(candidate, "non-Biological model with layers");

        KernelInterface iface;
        Atten_Mode original;
        original.absModel = OceanAbsorptionModel::Thorpe;
        iface.set_AttenUnit(original);
        test.requireThrows(
            [&] { iface.set_AttenUnit(candidate); },
            "invalid C++ setter input must throw");
        test.require(
            attenuationModesEqual(iface.getParams_const().AttenUnit, original),
            "failed C++ setter must preserve the original attenuation mode");
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

    void test_field_worker_honors_mode_limit(
        ook_test::TestRunner &test)
    {
        auto params = minimal_field_params();
        params.is_Velocity = false;
        params.SBP.isSet = false;
        params.MLimit = 1;
        params.Pos.Ro.setZero();

        auto eigen = minimal_eigen();
        auto limited_eigen = eigen;
        limited_eigen.M = 1;
        std::vector<std::complex<float>> expected(4, {0.0f, 0.0f});
        std::vector<std::complex<float>> unlimited(4, {0.0f, 0.0f});
        std::vector<std::complex<float>> actual(4, {0.0f, 0.0f});
        Evaluate(
            limited_eigen, params, 0, 0, expected.data(), nullptr, nullptr);
        Evaluate(eigen, params, 0, 0, unlimited.data(), nullptr, nullptr);

        bool second_mode_contributes = false;
        for (size_t index = 0; index < expected.size(); ++index)
        {
            if (std::abs(unlimited[index] - expected[index]) > 1.0e-5f)
            {
                second_mode_contributes = true;
            }
        }
        test.require(
            second_mode_contributes,
            "FieldWorker MLimit fixture must distinguish the second mode");

        OOK_output output{};
        output.eigen = &eigen;
        output.u_AllSources = actual.data();
        OpenOceanKraken::FieldWorker(0, params, output);

        for (size_t index = 0; index < actual.size(); ++index)
        {
            test.require(
                std::isfinite(actual[index].real()) &&
                    std::isfinite(actual[index].imag()),
                "FieldWorker MLimit pressure must be finite");
            test.requireNear(
                std::abs(actual[index] - expected[index]), 0.0, 2.0e-5,
                "FieldWorker MLimit=1 must sum only the first mode");
        }
        test.require(
            output.eigen[0].M == 2,
            "FieldWorker must not change the solved/exported mode count");
    }

    void test_field_worker_rejects_nonpositive_mode_limit(
        ook_test::TestRunner &test)
    {
        auto params = minimal_field_params();
        params.is_Velocity = false;
        params.SBP.isSet = false;
        params.MLimit = 0;
        auto eigen = minimal_eigen();
        std::vector<std::complex<float>> actual(4, {0.0f, 0.0f});
        OOK_output output{};
        output.eigen = &eigen;
        output.u_AllSources = actual.data();

        test.requireThrows(
            [&]() { OpenOceanKraken::FieldWorker(0, params, output); },
            "FieldWorker must preserve Evaluate's no-modes guard");
        test.require(
            output.eigen[0].M == 2,
            "rejected MLimit must not change the solved/exported mode count");
    }

    void test_evaluate_ad_identical_profiles_matches_single_profile(ook_test::TestRunner &test)
    {
        auto params = minimal_field_params();
        params.is_Velocity = false;
        params.SBP.isSet = false;
        params.Pos.Ro.setZero();
        params.NProf = 2;
        params.MLimit = 9999;
        params.RProf.resize(2);
        params.RProf << 0.0, 5.0;
        params.SSP.push_back(params.SSP.front());

        auto eigen = minimal_eigen();
        std::vector<EigenParams> profiles = {eigen, eigen};
        std::vector<std::complex<float>> expected(4, {0.0f, 0.0f});
        std::vector<std::complex<float>> actual(4, {0.0f, 0.0f});
        Evaluate(eigen, params, 0, 0, expected.data(), nullptr, nullptr);
        EvaluateAD(profiles.data(), static_cast<int>(profiles.size()), params, 0, actual.data());

        for (size_t i = 0; i < actual.size(); ++i)
        {
            test.requireNear(std::abs(actual[i] - expected[i]), 0.0, 2.0e-5,
                             "identical-profile EvaluateAD must match range-independent Evaluate");
        }
    }

    void test_evaluate_ad_modes_and_guards(ook_test::TestRunner &test)
    {
        auto params = minimal_field_params();
        params.is_Velocity = false;
        params.NProf = 2;
        params.MLimit = 9999;
        params.RProf.resize(2);
        params.RProf << 0.0, 5.0;
        params.SSP.push_back(params.SSP.front());
        auto eigen = minimal_eigen();
        std::vector<EigenParams> profiles = {eigen, eigen};
        std::vector<std::complex<float>> output(4, {0.0f, 0.0f});

        params.SBP.isSet = true;
        EvaluateAD(profiles.data(), 2, params, 0, output.data());
        test.require(std::isfinite(output.back().real()), "EvaluateAD beam-pattern output must be finite");

        auto bad_beam = params;
        bad_beam.SBP.theta.setZero();
        test.requireThrows([&]() { EvaluateAD(profiles.data(), 2, bad_beam, 0, output.data()); },
                           "EvaluateAD must reject duplicate beam-pattern angles");

        params.SBP.isSet = false;
        params.SourceType = Source_Mode::MODE_R_Point;
        params.coherenceType = CoherenceType::Incoherent;
        EvaluateAD(profiles.data(), 2, params, 0, output.data());
        test.require(std::abs(output[0]) == 0.0f, "EvaluateAD point source must vanish at zero range");

        auto three_profiles = params;
        three_profiles.coherenceType = CoherenceType::Coherent;
        three_profiles.NProf = 3;
        three_profiles.RProf.resize(3);
        three_profiles.RProf << 0.0, 3.0, 6.0;
        three_profiles.SSP.push_back(three_profiles.SSP.front());
        std::vector<EigenParams> eigen_three = {eigen, eigen, eigen};
        eigen_three[2].M = 1;
        EvaluateAD(eigen_three.data(), 3, three_profiles, 0, output.data());
        test.require(std::isfinite(output.back().real()), "EvaluateAD multi-segment output must be finite");

        test.requireThrows([&]() { EvaluateAD(nullptr, 2, params, 0, output.data()); },
                           "EvaluateAD must reject a null profile array");
        test.requireThrows([&]() { EvaluateAD(profiles.data(), 1, params, 0, output.data()); },
                           "EvaluateAD must reject a profile count mismatch");
        auto bad_ranges = params;
        bad_ranges.RProf(0) = 1.0;
        test.requireThrows([&]() { EvaluateAD(profiles.data(), 2, bad_ranges, 0, output.data()); },
                           "EvaluateAD must require RProf to start at zero");
        test.requireThrows([&]() { EvaluateAD(profiles.data(), 2, params, 1, output.data()); },
                           "EvaluateAD must reject an invalid source index");
        std::vector<EigenParams> incomplete = profiles;
        incomplete[1].M = 0;
        test.requireThrows([&]() { EvaluateAD(incomplete.data(), 2, params, 0, output.data()); },
                           "EvaluateAD must reject an incomplete eigen profile");
        auto decreasing_receivers = params;
        decreasing_receivers.Pos.Rr << 10.0, 5.0;
        test.requireThrows([&]() { EvaluateAD(profiles.data(), 2, decreasing_receivers, 0, output.data()); },
                           "EvaluateAD must reject decreasing receiver ranges");
        test.requireThrows([&]() { EvaluateAD(profiles.data(), 2, params, 0, nullptr); },
                           "EvaluateAD must reject a null pressure buffer");
        auto zero_limit = params;
        zero_limit.MLimit = 0;
        test.requireThrows([&]() { EvaluateAD(profiles.data(), 2, zero_limit, 0, output.data()); },
                           "EvaluateAD must reject an empty MLimit result");
    }

    void test_evaluate_cm_identical_profiles_preserve_amplitude(ook_test::TestRunner &test)
    {
        auto params = minimal_field_params();
        params.is_Velocity = false;
        params.SBP.isSet = false;
        params.SourceType = Source_Mode::MODE_X_Line;
        params.NProf = 2;
        params.MLimit = 1;
        params.RProf.resize(2);
        params.RProf << 0.0, 10.0;
        params.SSP.push_back(params.SSP.front());
        for (auto &ssp : params.SSP)
        {
            ssp.HSTop.BC = BC_Mode::MODE_V_Vacuum;
            ssp.HSBot.BC = BC_Mode::MODE_R_Rigid;
            ssp.rho.setOnes();
        }
        params.Pos.Rr << 0.0, 10.0;
        params.Pos.Rz << 0.0, 1.0;

        EigenParams eigen;
        eigen.M = 1;
        eigen.k.resize(1);
        eigen.k << std::complex<double>(0.2, 0.0);
        eigen.PsiS = Eigen::MatrixXcd::Ones(1, 1);
        eigen.PsiR = Eigen::MatrixXcd::Ones(1, 2);
        eigen.ModeZ.resize(2);
        eigen.ModeZ << 0.0, 1.0;
        eigen.PhiMode = Eigen::MatrixXcd::Ones(1, 2);
        std::vector<EigenParams> profiles = {eigen, eigen};

        std::vector<std::complex<float>> actual(4, {0.0f, 0.0f});
        EvaluateCM(profiles.data(), static_cast<int>(profiles.size()), params, 0, actual.data());

        constexpr float pi_single = 3.1415926f;
        const std::complex<float> i(0.0f, 1.0f);
        const std::complex<float> initial =
            std::sqrt(2.0f * pi_single) * std::exp(i * (pi_single / 4.0f)) / 0.2f;
        const std::complex<float> propagated =
            initial * std::exp(-i * std::complex<float>(0.2f, 0.0f) * 10.0f);
        for (int iz = 0; iz < 2; ++iz)
        {
            test.requireNear(std::abs(actual[GetFieldAddr(0, iz, 0, &params.Pos)] - initial), 0.0, 2.0e-5,
                             "EvaluateCM initial amplitude mismatch");
            test.requireNear(std::abs(actual[GetFieldAddr(0, iz, 1, &params.Pos)] - propagated), 0.0, 2.0e-5,
                             "EvaluateCM identical-profile crossing must preserve amplitude");
        }

        auto beamed = params;
        beamed.SBP.isSet = true;
        std::vector<std::complex<float>> beamed_output(4, {0.0f, 0.0f});
        EvaluateCM(profiles.data(), static_cast<int>(profiles.size()), beamed, 0, beamed_output.data());
        test.require(std::abs(beamed_output.back() - actual.back()) > 1.0e-5f,
                     "EvaluateCM source beam pattern must alter modal excitation");
    }

    void test_evaluate_cm_modes_and_guards(ook_test::TestRunner &test)
    {
        auto params = minimal_field_params();
        params.is_Velocity = false;
        params.SBP.isSet = false;
        params.SourceType = Source_Mode::MODE_R_Point;
        params.coherenceType = CoherenceType::Coherent;
        params.NProf = 2;
        params.MLimit = 1;
        params.RProf.resize(2);
        params.RProf << 0.0, 10.0;
        params.SSP.push_back(params.SSP.front());
        for (auto &ssp : params.SSP)
        {
            ssp.HSTop.BC = BC_Mode::MODE_V_Vacuum;
            ssp.HSBot.BC = BC_Mode::MODE_R_Rigid;
            ssp.rho.setOnes();
        }
        params.Pos.Rr << 0.0, 10.0;

        EigenParams eigen;
        eigen.M = 1;
        eigen.k.resize(1);
        eigen.k << std::complex<double>(0.2, 0.0);
        eigen.PsiS = Eigen::MatrixXcd::Ones(1, 1);
        eigen.PsiR = Eigen::MatrixXcd::Ones(1, 2);
        eigen.ModeZ.resize(2);
        eigen.ModeZ << 0.0, 1.0;
        eigen.PhiMode = Eigen::MatrixXcd::Ones(1, 2);
        std::vector<EigenParams> profiles = {eigen, eigen};
        std::vector<std::complex<float>> output(4, {0.0f, 0.0f});

        EvaluateCM(profiles.data(), 2, params, 0, output.data());
        test.require(std::isfinite(output.back().real()), "EvaluateCM point-source output must be finite");

        auto three_profiles = params;
        three_profiles.NProf = 3;
        three_profiles.RProf.resize(3);
        three_profiles.RProf << 0.0, 5.0, 10.0;
        three_profiles.SSP.push_back(three_profiles.SSP.front());
        std::vector<EigenParams> eigen_three = {eigen, eigen, eigen};
        EvaluateCM(eigen_three.data(), 3, three_profiles, 0, output.data());
        test.require(std::isfinite(output.back().real()), "EvaluateCM multi-segment output must be finite");

        test.requireThrows([&]() { EvaluateCM(nullptr, 2, params, 0, output.data()); },
                           "EvaluateCM must reject a null profile array");
        test.requireThrows([&]() { EvaluateCM(profiles.data(), 1, params, 0, output.data()); },
                           "EvaluateCM must reject a profile count mismatch");
        auto incoherent = params;
        incoherent.coherenceType = CoherenceType::Incoherent;
        test.requireThrows([&]() { EvaluateCM(profiles.data(), 2, incoherent, 0, output.data()); },
                           "EvaluateCM must reject incoherent propagation");
        test.requireThrows([&]() { EvaluateCM(profiles.data(), 2, params, 1, output.data()); },
                           "EvaluateCM must reject an invalid source index");
        std::vector<EigenParams> incomplete = profiles;
        incomplete[1].PhiMode.resize(0, 0);
        test.requireThrows([&]() { EvaluateCM(incomplete.data(), 2, params, 0, output.data()); },
                           "EvaluateCM must reject an incomplete eigen profile");
        test.requireThrows([&]() { EvaluateCM(profiles.data(), 2, params, 0, nullptr); },
                           "EvaluateCM must reject a null pressure buffer");
        auto zero_limit = params;
        zero_limit.MLimit = 0;
        test.requireThrows([&]() { EvaluateCM(profiles.data(), 2, zero_limit, 0, output.data()); },
                           "EvaluateCM must reject an empty MLimit result");
        auto bad_profile_ranges = params;
        bad_profile_ranges.RProf << 0.0, 0.0;
        test.requireThrows([&]() { EvaluateCM(profiles.data(), 2, bad_profile_ranges, 0, output.data()); },
                           "EvaluateCM must reject non-increasing profile ranges");
        auto decreasing_receivers = params;
        decreasing_receivers.Pos.Rr << 10.0, 5.0;
        test.requireThrows([&]() { EvaluateCM(profiles.data(), 2, decreasing_receivers, 0, output.data()); },
                           "EvaluateCM must reject decreasing receiver ranges");
        auto missing_ssp = params;
        missing_ssp.SSP.pop_back();
        test.requireThrows([&]() { EvaluateCM(profiles.data(), 2, missing_ssp, 0, output.data()); },
                           "EvaluateCM must require one SSP per profile");
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

        Atten_Mode biological;
        biological.absModel = OceanAbsorptionModel::Biological;
        biological.biologicalLayers = {{0.0, std::numeric_limits<double>::max(), 1000.0, 5.0, 0.04}};
        HSInfo bio_top = top;
        HSInfo bio_bottom = bot;
        HSInfo baseline_top = top;
        HSInfo baseline_bottom = bot;
        Atten_Mode no_volume;
        OpenOceanKraken::ssp::UpdateHSLoss(freq, freq0, biological, bio_top, bio_bottom);
        OpenOceanKraken::ssp::UpdateHSLoss(freq, freq0, no_volume, baseline_top, baseline_bottom);
        test.requireNear(std::imag(bio_top.cp), std::imag(baseline_top.cp), 1.0e-15, "top half-space must exclude Biological attenuation");
        test.requireNear(std::imag(bio_bottom.cp), std::imag(baseline_bottom.cp), 1.0e-15, "bottom half-space must exclude Biological attenuation");

        auto bio_ssp = ssp_fixture();
        auto baseline_ssp = ssp_fixture();
        const Eigen::Index stored_points = bio_ssp.z.size();
        const std::complex<double> nan(std::numeric_limits<double>::quiet_NaN(), 0.0);
        bio_ssp.cpCoef.setConstant(nan);
        bio_ssp.csCoef.setConstant(nan);
        bio_ssp.cpSpline.setConstant(nan);
        bio_ssp.csSpline.setConstant(nan);
        Atten_Mode whole_water_bio;
        whole_water_bio.absModel = OceanAbsorptionModel::Biological;
        whole_water_bio.biologicalLayers = {{0.0, 100.0, 50.0, 5.0, 0.04}};
        Atten_Mode no_volume_ssp;
        OpenOceanKraken::ssp::UpdateSSPLoss(50.0, 50.0, baseline_ssp.NMedia, SSP_Mode::MODE_P_cPCHIP, no_volume_ssp, baseline_ssp);
        OpenOceanKraken::ssp::UpdateSSPLoss(50.0, 50.0, bio_ssp.NMedia, SSP_Mode::MODE_P_cPCHIP, whole_water_bio, bio_ssp);
        test.require(bio_ssp.z.size() == stored_points, "Biological attenuation must not insert SSP boundary nodes");
        test.require(std::imag(bio_ssp.cp(1)) > std::imag(baseline_ssp.cp(1)) && std::imag(bio_ssp.cs(1)) > std::imag(baseline_ssp.cs(1)), "P and non-zero S speed must use the same Biological mode");
        test.require(bio_ssp.cpCoef.allFinite() && bio_ssp.csCoef.allFinite(), "PCHIP coefficients must be rebuilt after Biological loss");
        OpenOceanKraken::ssp::UpdateSSPLoss(50.0, 50.0, bio_ssp.NMedia, SSP_Mode::MODE_S_cCubic, whole_water_bio, bio_ssp);
        test.require(bio_ssp.cpSpline.allFinite() && bio_ssp.csSpline.allFinite(), "spline coefficients must be rebuilt after Biological loss");

        OOK_parameters params = minimal_field_params();
        EigenParams eigen = minimal_eigen();
        std::vector<std::complex<float>> u(4), v(4), h(4);
        Evaluate(eigen, params, 0, 0, u.data(), v.data(), h.data());
        test.require(std::isfinite(std::real(u[0])), "coherent field output must be finite");

        EigenParams directional_eigen = minimal_eigen();
        directional_eigen.dPsidzR.setZero();
        std::fill(u.begin(), u.end(), std::complex<float>{});
        std::fill(v.begin(), v.end(), std::complex<float>{});
        std::fill(h.begin(), h.end(), std::complex<float>{});
        Evaluate(directional_eigen, params, 0, 0, u.data(), v.data(), h.data());
        test.require(std::abs(v[0]) < 1.0e-7f, "coherent vertical velocity must use dPsiR/dz");
        test.require(std::abs(h[0]) > 1.0e-7f, "coherent horizontal velocity must use PsiR");

        params.SourceType = Source_Mode::MODE_R_Point;
        params.coherenceType = CoherenceType::Incoherent;
        std::fill(u.begin(), u.end(), std::complex<float>{});
        std::fill(v.begin(), v.end(), std::complex<float>{});
        std::fill(h.begin(), h.end(), std::complex<float>{});
        Evaluate(eigen, params, 0, 0, u.data(), v.data(), h.data());
        test.require(std::isfinite(std::real(u[1])), "incoherent field output must be finite");

        params.SourceType = Source_Mode::MODE_X_Line;
        std::fill(u.begin(), u.end(), std::complex<float>{});
        std::fill(v.begin(), v.end(), std::complex<float>{});
        std::fill(h.begin(), h.end(), std::complex<float>{});
        Evaluate(directional_eigen, params, 0, 0, u.data(), v.data(), h.data());
        test.require(std::abs(v[0]) < 1.0e-7f, "incoherent vertical velocity must use dPsiR/dz");
        test.require(std::abs(h[0]) > 1.0e-7f, "incoherent horizontal velocity must use PsiR");

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

        paramsBase *input_base = new input_Freq();
        delete input_base;
        input_base = new input_Sz_Rz_RR();
        delete input_base;
        input_base = new input_reflcoef();
        delete input_base;
        input_base = new input_sbp();
        delete input_base;
        input_base = new input_SSP();
        delete input_base;

        outputBase *output_base = new OpenOceanKraken::output_Eigen();
        delete output_base;
        output_base = new OpenOceanKraken::output_Field();
        delete output_base;
        output_base = new outputBase();
        delete output_base;

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
        OpenOceanKraken::ssp::SSPLayer copied_layer;
        copied_layer = layer;
        OpenOceanKraken::ssp::SSPLayer moved_layer;
        moved_layer = std::move(copied_layer);
        test.require(moved_layer.npts == layer.npts, "SSPLayer move assignment point count mismatch");
        test.requireNear(moved_layer.z[moved_layer.npts - 1], layer.z[layer.npts - 1], 1.0e-12, "SSPLayer move assignment depth mismatch");
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

        Eigen::MatrixXd matrix(2, 2);
        matrix << 1.0, 2.0, 3.0, 4.0;
        OpenOcean_json matrix_json = matrix;
        const auto decoded_matrix = matrix_json.get<Eigen::MatrixXd>();
        test.requireNear(decoded_matrix(1, 0), 3.0, 1.0e-12, "2D Eigen JSON roundtrip mismatch");

        OpenOcean_json ragged_json = OpenOcean_json::parse("[[1.0, 2.0], [3.0]]");
        test.requireThrows([&]() { (void)ragged_json.get<Eigen::MatrixXd>(); }, "ragged Eigen JSON matrix must throw");
        OpenOcean_json wrong_columns = OpenOcean_json::parse("[[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]");
        test.requireThrows([&]() { (void)wrong_columns.get<Eigen::Matrix<double, Eigen::Dynamic, 2>>(); }, "fixed-column Eigen JSON mismatch must throw");

        Eigen::VectorXd merge_x(3), merge_y(1), merged;
        Eigen::VectorXi ix, iy;
        int merged_count = 0;
        merge_x << 1.0, 2.0, 3.0;
        merge_y << 1.0;
        OpenOceanKraken::MergeVectors(merge_x, merge_y, merged, merged_count, ix, iy);
        test.require(merged_count == 3 && ix[2] == 2 && iy[0] == 0, "MergeVectors exhausted-y path mismatch");

        merge_x.resize(1);
        merge_y.resize(1);
        merge_x << 1.0 + std::numeric_limits<double>::epsilon();
        merge_y << 1.0;
        OpenOceanKraken::MergeVectors(merge_x, merge_y, merged, merged_count, ix, iy);
        test.require(merged_count == 1 && ix[0] == 0 && iy[0] == 0, "MergeVectors near-duplicate x path mismatch");

        Eigen::VectorXd one(1), tab(3), weights(3);
        Eigen::VectorXi weight_indices(3);
        one << 5.0;
        tab << 4.0, 5.0, 6.0;
        weights.setConstant(123.0);
        weight_indices.setConstant(-7);
        OpenOceanKraken::Weight_dble(one, 1, tab, 3, weights, weight_indices);
        for (int query = 0; query < 3; ++query)
        {
            test.requireNear(
                weights[query], 0.0, 1.0e-12,
                "Weight_dble singleton weight mismatch at query " +
                    std::to_string(query));
            test.require(
                weight_indices[query] == 0,
                "Weight_dble singleton index mismatch at query " +
                    std::to_string(query));
        }
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

            KernelInterface iface;
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

            KernelInterface iface;
            const auto env_path = entry.path();
            const auto name = env_path.filename().string();
            if (name == "option_default_chars.env")
            {
                test.require(!iface.from_env(env_path.string()),
                             name + ": unknown fourth TOP OPTION character must fail");
                ++parsed;
                continue;
            }
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
                ThreadPool pool(1);
                iface.setThreadPool(pool);
                iface.setNumThreads(1);
                iface.runEigen();
                test.require(iface.getOutput_const().eigen[0].M > 0, name + ": acoustic top half-space must produce modes");
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
            KernelInterface json_iface;
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
            KernelInterface brc_iface;
            test.require(!brc_iface.from_env(brc_env.string()), "KRAKEN bottom F/.brc must fail fast");
        }

        const fs::path irc_env = test_root / "neggradK_irc.env";
        const fs::path irc_flp = test_root / "neggradK_irc.flp";
        if (!fs::exists(irc_env) || !fs::exists(irc_flp))
        {
            return;
        }

        KernelInterface missing_irc_iface;
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

        KernelInterface irc_iface;
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
        KernelInterface iface;
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
        KernelInterface env_iface;
        test.require(env_iface.from_env((test_root / "MunkK.env").string()), "MunkK from_env before JSON must succeed");
        const fs::path json_path = fs::temp_directory_path() / "ook_munk_roundtrip.json";
        test.require(env_iface.to_json(json_path.string()), "to_json must succeed for MunkK");

        KernelInterface json_iface;
        test.require(json_iface.from_json(json_path.string()), "from_json must load MunkK snapshot");
        const auto &params = json_iface.getParams_const();
        test.requireNear(params.freqinfo.freq, 50.0, 1.0e-9, "JSON MunkK freq mismatch");
        test.require(params.SSP.size() == params.sspInput.size(), "JSON MunkK SSP conversion count mismatch");
        test.require(params.Pos.NRz == 501, "JSON MunkK NRz mismatch");
        fs::remove(json_path);
    }

    void test_json_roundtrip_from_params(ook_test::TestRunner &test)
    {
        KernelInterface params_iface;
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

        KernelInterface json_iface;
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
        test.require(
            !base_json["AttenUnit"].contains("BiologicalLayers"),
            "non-Biological JSON output must omit BiologicalLayers");

        Atten_Mode biological;
        biological.absModel = OceanAbsorptionModel::Biological;
        biological.biologicalLayers = {
            {10.0, 30.0, 1000.0, 5.0, 0.04},
            {40.0, 60.0, 1200.0, 4.0, 0.02}};
        params_iface.set_AttenUnit(biological);
        const auto biological_json =
            OpenOcean_json::parse(params_iface.to_json_string());
        test.require(
            biological_json["AttenUnit"]["OceanAbsorptionModel"] ==
                "Biological",
            "JSON must emit Biological model");
        test.require(
            biological_json["AttenUnit"]["BiologicalLayers"].size() == 2,
            "JSON must preserve Biological layer count");
        test.requireNear(
            biological_json["AttenUnit"]["BiologicalLayers"][1]["f0"]
                .get<double>(),
            1200.0, 0.0, "JSON must preserve layer order and fields");

        const fs::path biological_params_json =
            fs::temp_directory_path() / "ook_biological_params_roundtrip.json";
        test.require(
            params_iface.to_json(biological_params_json.string()),
            "Biological params must export to JSON");
        KernelInterface biological_params_copy;
        test.require(
            biological_params_copy.from_json(biological_params_json.string()),
            "Biological params JSON must reload");
        test.require(
            attenuationModesEqual(
                params_iface.getParams_const().AttenUnit,
                biological_params_copy.getParams_const().AttenUnit),
            "JSON must preserve every Biological layer field and order");
        fs::remove(biological_params_json);

        Atten_Mode empty_biological = biological;
        empty_biological.biologicalLayers.clear();
        params_iface.set_AttenUnit(empty_biological);
        const auto empty_biological_json =
            OpenOcean_json::parse(params_iface.to_json_string());
        test.require(
            empty_biological_json["AttenUnit"]["BiologicalLayers"].empty(),
            "Biological JSON must preserve an empty layer array");
        const fs::path empty_biological_path =
            fs::temp_directory_path() / "ook_empty_biological_roundtrip.json";
        write_text(empty_biological_path, empty_biological_json.dump(2));
        KernelInterface empty_biological_copy;
        test.require(
            empty_biological_copy.from_json(empty_biological_path.string()),
            "empty Biological layer array must reload");
        test.require(
            attenuationModesEqual(
                empty_biological,
                empty_biological_copy.getParams_const().AttenUnit),
            "empty Biological layer array must round-trip");
        fs::remove(empty_biological_path);

        const fs::path bio_roundtrip_root =
            fs::temp_directory_path() / "ook_biological_env_json_roundtrip";
        fs::create_directories(bio_roundtrip_root);
        const fs::path valid_bio_env = bio_roundtrip_root / "case.env";
        write_text(
            valid_bio_env,
            biological_env_text(
                "CVWB", "1\n20 40 1000 5 0.04\n"));
        write_text(bio_roundtrip_root / "case.flp", single_profile_flp);
        KernelInterface env_source;
        test.require(
            env_source.from_env(valid_bio_env.string()),
            "Biological ENV must load before JSON round-trip");
        const fs::path env_json =
            fs::temp_directory_path() / "ook_biological_env_roundtrip.json";
        test.require(
            env_source.to_json(env_json.string()),
            "Biological ENV must export to JSON");
        KernelInterface env_json_copy;
        test.require(
            env_json_copy.from_json(env_json.string()),
            "Biological ENV-derived JSON must reload");
        test.require(
            attenuationModesEqual(
                env_source.getParams_const().AttenUnit,
                env_json_copy.getParams_const().AttenUnit),
            "ENV to JSON must preserve every layer field and order");
        fs::remove(env_json);
        fs::remove_all(bio_roundtrip_root);

        auto expect_bad_json = [&](OpenOcean_json bad, const std::string &suffix) {
            const fs::path bad_path = fs::temp_directory_path() / ("ook_bad_json_" + suffix + ".json");
            {
                std::ofstream bad_file(bad_path);
                bad_file << bad.dump(2);
            }
            KernelInterface bad_iface;
            Atten_Mode thorpe;
            thorpe.absModel = OceanAbsorptionModel::Thorpe;
            bad_iface.set_AttenUnit(thorpe);
            test.require(!bad_iface.from_json(bad_path.string()), "bad JSON variant must fail: " + suffix);
            test.require(
                attenuationModesEqual(
                    bad_iface.getParams_const().AttenUnit, thorpe),
                "failed JSON import must preserve Thorpe attenuation mode: " + suffix);
            fs::remove(bad_path);
        };

        const fs::path non_biological_empty_path =
            fs::temp_directory_path() / "ook_non_biological_empty_layers.json";
        auto non_biological_empty = base_json;
        non_biological_empty["AttenUnit"]["BiologicalLayers"] =
            OpenOcean_json::array();
        write_text(non_biological_empty_path, non_biological_empty.dump(2));
        KernelInterface non_biological_empty_copy;
        test.require(
            non_biological_empty_copy.from_json(non_biological_empty_path.string()),
            "non-Biological JSON with empty BiologicalLayers must reload");
        const auto normalized_non_biological_json = OpenOcean_json::parse(
            non_biological_empty_copy.to_json_string());
        test.require(
            !normalized_non_biological_json["AttenUnit"].contains(
                "BiologicalLayers"),
            "non-Biological JSON must normalize away empty BiologicalLayers");
        fs::remove(non_biological_empty_path);

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

        bad = biological_json;
        bad["AttenUnit"].erase("BiologicalLayers");
        expect_bad_json(bad, "biological_layers_missing");
        bad = biological_json;
        bad["AttenUnit"]["BiologicalLayers"] = 12;
        expect_bad_json(bad, "biological_layers_type");
        bad = biological_json;
        bad["AttenUnit"]["BiologicalLayers"][0].erase("f0");
        expect_bad_json(bad, "biological_layer_field_missing");
        bad = biological_json;
        bad["AttenUnit"]["BiologicalLayers"][0]["f0"] = "not-a-number";
        expect_bad_json(bad, "biological_layer_field_type");
        bad = biological_json;
        bad["AttenUnit"]["BiologicalLayers"] = OpenOcean_json::array();
        for (int i = 0; i <= MaxBioLayers; ++i)
        {
            bad["AttenUnit"]["BiologicalLayers"].push_back(
                {{"Z1", 10.0}, {"Z2", 30.0}, {"f0", 1000.0},
                 {"Q", 5.0}, {"a0", 0.04}});
        }
        expect_bad_json(bad, "biological_layers_over_limit");
        bad = base_json;
        bad["AttenUnit"]["BiologicalLayers"] =
            biological_json["AttenUnit"]["BiologicalLayers"];
        expect_bad_json(bad, "non_biological_layers_nonempty");
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
                KernelInterface iface;
                test.require(!iface.from_env(path.string()), name + ": bad ENV fixture must fail from_env");
                ++env_cases;
            }
            else if (path.extension() == ".json" && name != "_manifest.json")
            {
                KernelInterface iface;
                test.require(!iface.from_json(path.string()), name + ": bad JSON fixture must fail from_json");
                ++json_cases;
            }
        }

        test.require(env_cases > 0, "LCOV bad input coverage requires at least one bad ENV fixture");
        test.require(json_cases > 0, "LCOV bad input coverage requires at least one bad JSON fixture");
    }

    void test_mult_profile_velocity_is_rejected_until_supported(ook_test::TestRunner &test, const fs::path &test_root)
    {
        fs::path env_path = test_root / "solve3_mode_loss.env";
        if (!fs::exists(env_path))
        {
            env_path = test_root / "solve3_examples" / "solve3_mode_loss.env";
        }
        if (!fs::exists(env_path))
        {
            return;
        }

        ThreadPool pool(1);
        KernelInterface iface(pool);
        iface.setNumThreads(1);
        test.require(iface.from_env(env_path.string()), "multi-profile field guard fixture must load");
        test.require(iface.getParams_const().SSP.size() > 1, "multi-profile field guard fixture must contain multiple SSP profiles");
        iface.getParams().is_Velocity = true;
        iface.runEigen();

        bool rejected = false;
        try
        {
            iface.runField();
        }
        catch (const std::logic_error &error)
        {
            rejected = std::string(error.what()) == "multi-profile velocity is not implemented";
        }
        test.require(rejected, "runField must explicitly reject unsupported multi-profile velocity");
    }

    void test_mult_profile_coupled_field_runs(ook_test::TestRunner &test, const fs::path &test_root)
    {
        for (const std::string &case_name : {"solve3_mode_loss", "solve3_mode_gain"})
        {
            const fs::path env_path = test_root / (case_name + ".env");
            if (!fs::exists(env_path))
            {
                continue;
            }

            ThreadPool pool(1);
            KernelInterface iface(pool);
            iface.setNumThreads(1);
            test.require(iface.from_env(env_path.string()), case_name + ": multi-profile coupled fixture must load");
            test.require(iface.getParams_const().modeType == ModeType::Couple, case_name + ": coupled FLP mode mismatch");
            iface.runEigen();
            if (case_name == "solve3_mode_gain")
            {
                const auto &output = iface.getOutput_const();
                test.require(output.eigen[1].M > output.eigen[0].M,
                             "solve3_mode_gain must increase the available mode count");
            }
            iface.runField();

            const auto &params = iface.getParams_const();
            const std::complex<float> *pressure = iface.get_u_AllSources();
            const size_t count = static_cast<size_t>(params.Pos.NSz) *
                                 static_cast<size_t>(params.Pos.NRz_per_range) *
                                 static_cast<size_t>(params.Pos.NRr);
            for (size_t i = 0; i < count; ++i)
            {
                test.require(std::isfinite(pressure[i].real()) && std::isfinite(pressure[i].imag()),
                             case_name + ": multi-profile coupled pressure must be finite");
            }
        }
    }

    void test_mult_profile_adiabatic_field_runs(ook_test::TestRunner &test, const fs::path &test_root)
    {
        const fs::path env_path = test_root / "adiabatic_mode_loss.env";
        if (!fs::exists(env_path))
        {
            return;
        }

        ThreadPool pool(1);
        KernelInterface iface(pool);
        iface.setNumThreads(1);
        test.require(iface.from_env(env_path.string()), "multi-profile adiabatic fixture must load");
        test.require(iface.getParams_const().modeType == ModeType::Adiabatic, "adiabatic FLP mode mismatch");
        iface.runEigen();
        iface.runField();

        const auto &params = iface.getParams_const();
        const std::complex<float> *pressure = iface.get_u_AllSources();
        const size_t count = static_cast<size_t>(params.Pos.NSz) *
                             static_cast<size_t>(params.Pos.NRz_per_range) *
                             static_cast<size_t>(params.Pos.NRr);
        for (size_t i = 0; i < count; ++i)
        {
            test.require(std::isfinite(pressure[i].real()) && std::isfinite(pressure[i].imag()),
                         "multi-profile adiabatic pressure must be finite");
        }
    }

    void test_mult_profile_flp_metadata_is_serialized(ook_test::TestRunner &test, const fs::path &test_root)
    {
        const fs::path env_path = test_root / "solve3_mode_loss.env";
        if (!fs::exists(env_path))
        {
            return;
        }

        KernelInterface iface;
        test.require(iface.from_env(env_path.string()), "multi-profile FLP metadata fixture must load");
        const auto &params = iface.getParams_const();
        test.require(params.MLimit == 9999, "multi-profile FLP MLimit mismatch");
        test.require(params.NProf == 2, "multi-profile FLP NProf mismatch");
        test.require(params.RProf.size() == 2, "multi-profile FLP RProf size mismatch");
        test.requireNear(params.RProf(0), 0.0, 1.0e-9, "multi-profile first range mismatch");
        test.requireNear(params.RProf(1), 5000.0, 1.0e-9, "multi-profile second range mismatch");
        const std::string json = iface.to_json_string();
        test.require(json.find("\"MLimit\"") != std::string::npos, "multi-profile JSON must preserve MLimit");
        test.require(json.find("\"NProf\"") != std::string::npos, "multi-profile JSON must preserve NProf");
        test.require(json.find("\"RProf\"") != std::string::npos, "multi-profile JSON must preserve RProf");

        const fs::path json_path = fs::temp_directory_path() / "ook_mult_profile_roundtrip.json";
        test.require(iface.to_json(json_path.string()), "multi-profile JSON export must succeed");
        KernelInterface roundtrip;
        test.require(roundtrip.from_json(json_path.string()), "multi-profile JSON import must succeed");
        const auto &roundtrip_params = roundtrip.getParams_const();
        test.require(roundtrip_params.MLimit == 9999, "roundtrip MLimit mismatch");
        test.require(roundtrip_params.NProf == 2, "roundtrip NProf mismatch");
        test.requireNear(roundtrip_params.RProf(1), 5000.0, 1.0e-9, "roundtrip second profile range mismatch");
        fs::remove(json_path);
    }

    void test_invalid_mult_profile_flp_is_rejected(ook_test::TestRunner &test, const fs::path &test_root)
    {
        const fs::path source_env = test_root / "solve3_mode_loss.env";
        if (!fs::exists(source_env))
        {
            return;
        }

        const fs::path root = fs::temp_directory_path() / "ook_invalid_mult_profile_flp";
        fs::create_directories(root);
        const fs::path env_path = root / "case.env";
        const fs::path flp_path = root / "case.flp";
        fs::copy_file(source_env, env_path, fs::copy_options::overwrite_existing);

        auto write_flp = [&flp_path](int nprof, const std::string &rprof)
        {
            std::ofstream flp(flp_path);
            flp << "/,\n'RC'\n9999,\n" << nprof << "\n" << rprof
                << "\n101\n0.0 5.0 /\n1\n50.0 /\n101\n0.0 100.0 /\n101\n0.0 0.0 /\n";
        };

        write_flp(2, "0.0 0.0 /");
        KernelInterface non_increasing;
        test.require(!non_increasing.from_env(env_path.string()), "non-increasing RProf must fail from_env");

        write_flp(3, "0.0 5.0 10.0 /");
        KernelInterface count_mismatch;
        test.require(!count_mismatch.from_env(env_path.string()), "NProf/ENV profile mismatch must fail from_env");
        fs::remove_all(root);
    }

    void test_normalize_file_top_halfspace_bottom_wavenumbers(ook_test::TestRunner &test, const fs::path &test_root)
    {
        const fs::path env_path = test_root / "normalize_boundary_top_file.env";
        if (!fs::exists(env_path))
        {
            return;
        }

        ThreadPool pool(1);
        KernelInterface iface(pool);
        iface.setNumThreads(1);
        test.require(iface.from_env(env_path.string()), "Normalize F/A regression fixture must load");
        iface.runEigen();

        const auto &eigen = iface.getOutput_const().eigen[0];
        test.require(eigen.M == 7, "Normalize F/A regression must retain seven modes");
        const std::complex<double> kraken_k[] = {
            {0.20901788771152496, 0.00043763950816355646},
            {0.20555116236209870, 0.00135090237017720940},
            {0.19825963675975800, 0.00232429802417755130},
            {0.18705640733242035, 0.00269085052423179150},
            {0.16757051646709442, 0.00227311137132346630},
            {0.13855974376201630, 0.00241931318305432800},
            {0.09136671572923660, 0.00316715380176901800},
        };
        for (int mode = 0; mode < eigen.M; ++mode)
        {
            test.requireNear(std::real(eigen.k(mode)), std::real(kraken_k[mode]), 1.0e-7,
                             "Normalize F/A real wavenumber mismatch at mode " + std::to_string(mode + 1));
            test.requireNear(std::imag(eigen.k(mode)), std::imag(kraken_k[mode]), 1.0e-7,
                             "Normalize F/A imaginary wavenumber mismatch at mode " + std::to_string(mode + 1));
        }
    }

    void test_multilayer_elastic_stack_terminal_mode(ook_test::TestRunner &test, const fs::path &test_root)
    {
        ThreadPool pool(1);
        KernelInterface iface(pool);
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
        KernelInterface iface(pool);
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
        KernelInterface iface(pool);
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
        test_interface_velocity_and_export_workflow(test, env_root);
        test_top_bottom_line(test);
        test_interpolation_helpers(test);
        test_attenuation_and_reflection_helpers(test);
        test_biological_attenuation_validation(test);
        test_boundary_and_scatter_helpers(test);
        test_field_worker_rejects_nonpositive_mode_limit(test);
        test_field_worker_honors_mode_limit(test);
        test_evaluate_ad_identical_profiles_matches_single_profile(test);
        test_evaluate_ad_modes_and_guards(test);
        test_evaluate_cm_identical_profiles_preserve_amplitude(test);
        test_evaluate_cm_modes_and_guards(test);
        test_ssp_field_and_file_helpers(test);
        test_input_module_helpers(test);
        test_util_helpers(test);
        test_singleton_mode_position_interpolation(test);
        test_biological_env_transport(test);
        test_all_env_flp_cases(test, env_root);
        test_option_fixture_cases(test, option_root);
        test_kraken_reflection_option_guards(test, env_root);
        test_munk_exact_values(test, env_root);
        test_json_roundtrip_from_env(test, env_root);
        test_json_roundtrip_from_params(test);
        test_bad_input_cases(test, bad_root);
        test_mult_profile_flp_metadata_is_serialized(test, test_root);
        test_invalid_mult_profile_flp_is_rejected(test, test_root);
        test_mult_profile_velocity_is_rejected_until_supported(test, test_root);
        test_mult_profile_coupled_field_runs(test, test_root);
        test_mult_profile_adiabatic_field_runs(test, test_root);
        test_normalize_file_top_halfspace_bottom_wavenumbers(test, test_root);
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
