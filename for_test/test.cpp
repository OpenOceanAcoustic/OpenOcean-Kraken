#include "test.h"
#include "OpenOceanKrakenInterface.h"

#include <Eigen/Dense>
#include <algorithm>
#include <filesystem>
#include <fstream>
#include <vector>

namespace fs = std::filesystem;
using OpenOceanKraken::BC_Mode;
using OpenOceanKraken::Grid_Mode;
using OpenOceanKraken::Interface;
using OpenOceanKraken::OOK_parameters;
using OpenOceanKraken::Run_Mode;
using OpenOceanKraken::ssp::Range_Independent_Area;

namespace
{
    std::vector<std::string> expected_env_cases()
    {
        return {};
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
        params_iface.set_Sz(12.0, 24.0, 3);
        params_iface.set_Rz(0.0, 90.0, 4);
        params_iface.set_Rr(0.0, 3000.0, 4);
        params_iface.set_RunMode(Run_Mode::MODE_M_Modes);

        const fs::path json_path = fs::temp_directory_path() / "ook_params_roundtrip.json";
        test.require(params_iface.to_json(json_path.string()), "params to_json must succeed");

        Interface json_iface;
        test.require(json_iface.from_json(json_path.string()), "params from_json must succeed");
        const auto &params = json_iface.getParams_const();
        test.require(params.Title == title, "JSON params title mismatch");
        test.requireNear(params.freqinfo.freq, 77.0, 1.0e-9, "JSON params freq mismatch");
        test.require(params.Pos.NSz == 3, "JSON params NSz mismatch");
        test.require(params.Pos.NRz == 4, "JSON params NRz mismatch");
        test.require(params.Pos.NRr == 4, "JSON params NRr mismatch");
        test.require(params.runMode == Run_Mode::MODE_M_Modes, "JSON params run mode mismatch");
        test.require(params.SSP.size() == params.sspInput.size(), "JSON params SSP conversion count mismatch");
        fs::remove(json_path);
    }
}

int main(int argc, char **argv)
{
    ook_test::TestRunner test;
    try
    {
        const fs::path test_root = (argc > 1) ? fs::path(argv[1]) : fs::path("test");
        test.require(fs::exists(test_root), "test root must exist");
        test_lifecycle_and_guards(test);
        test_top_bottom_line(test);
        test_all_env_flp_cases(test, test_root);
        test_munk_exact_values(test, test_root);
        test_json_roundtrip_from_env(test, test_root);
        test_json_roundtrip_from_params(test);
        std::cout << "OOK interface tests passed. checks=" << test.checks() << std::endl;
        std::cout << "ENV/FLP coverage=all discovered ENV cases, JSON coverage=2/2" << std::endl;
        return 0;
    }
    catch (const std::exception &e)
    {
        std::cerr << "OOK interface tests failed: " << e.what() << std::endl;
        std::cerr << "checks before failure=" << test.checks() << std::endl;
        return 1;
    }
}
