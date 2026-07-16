#include "OpenOceanKrakencInterface.h"
#include "module/json_in_out.hpp"

#include <cassert>
#include <filesystem>
#include <fstream>
#include <stdexcept>

using namespace OpenOceanKrakenc;

#undef assert
#define assert(condition) do { if (!(condition)) throw std::runtime_error("assertion failed: " #condition); } while (false)

namespace
{
template <typename Enum>
void roundTripEnum(Enum value)
{
    OpenOcean_json encoded;
    to_json(encoded, value);
    Enum restored{};
    from_json(encoded, restored);
    assert(restored == value);
}
}

int main()
{
    for (SSP_Mode value : {SSP_Mode::MODE_N_n2Linear, SSP_Mode::MODE_C_cLinear,
                           SSP_Mode::MODE_P_cPCHIP, SSP_Mode::MODE_S_cCubic,
                           SSP_Mode::MODE_A_Analytic}) roundTripEnum(value);
    for (AttenuationUnit value : {AttenuationUnit::MODE_F_dB_per_m_kHz,
                                  AttenuationUnit::MODE_L_params_lose,
                                  AttenuationUnit::MODE_M_dB_per_m,
                                  AttenuationUnit::MODE_m_dB_per_m,
                                  AttenuationUnit::MODE_N_Nepers_per_m,
                                  AttenuationUnit::MODE_Q_Quality_Factor,
                                  AttenuationUnit::MODE_W_db_per_lambda}) roundTripEnum(value);
    for (BC_Mode value : {BC_Mode::MODE_R_Rigid, BC_Mode::MODE_V_Vacuum,
                          BC_Mode::MODE_F_File, BC_Mode::MODE_A_Half_space,
                          BC_Mode::MODE_G_Grain, BC_Mode::MODE_P_Precomputed}) roundTripEnum(value);
    for (Source_Mode value : {Source_Mode::MODE_R_Point, Source_Mode::MODE_X_Line,
                              Source_Mode::MODE_S_ScaledCylindrical}) roundTripEnum(value);
    roundTripEnum(Grid_Mode::MODE_R_Rectangular);
    roundTripEnum(Grid_Mode::MODE_I_Irregular);
    roundTripEnum(Run_Mode::MODE_M_Modes);
    roundTripEnum(Run_Mode::MODE_F_Field);
    roundTripEnum(Run_Mode::MODE_B_Both);
    roundTripEnum(CoherenceType::Coherent);
    roundTripEnum(CoherenceType::Incoherent);
    roundTripEnum(ModeType::Adiabatic);
    roundTripEnum(ModeType::Couple);
    roundTripEnum(OceanAbsorptionModel::None);
    roundTripEnum(OceanAbsorptionModel::Thorpe);
    roundTripEnum(OceanAbsorptionModel::FrancGarr);
    roundTripEnum(Media_Mode::MODE_A_Acoustic);
    roundTripEnum(Media_Mode::MODE_E_Elastic);

    const std::filesystem::path workspace = OPENOCEANKRAKENC_WORKSPACE_DIR;
    const std::filesystem::path env = workspace / "test" / "MunkKleaky.env";
    const std::filesystem::path output =
        std::filesystem::temp_directory_path() / "openocean_krakenc_roundtrip.json";
    const std::filesystem::path invalid =
        std::filesystem::temp_directory_path() / "openocean_krakenc_invalid.json";

    Interface source;
    assert(source.from_env(env.string()));
    source.set_Velocity_enable(true);
    const std::string encoded = source.to_json_string();
    const OpenOcean_json document = OpenOcean_json::parse(encoded);
    assert(document.at("schema") == "OpenOcean-Krakenc.parameters");
    assert(document.at("schema_version") == 1);
    assert(document.contains("Title"));
    assert(document.contains("freqinfo"));
    assert(document.contains("Pos"));
    assert(document.contains("sspInput"));
    assert(document.contains("ReflectionCoef"));
    assert(document.contains("SBP"));
    assert(document.at("AttenUnit").contains("AttenuationUnit"));
    assert(document.at("AttenUnit").contains("OceanAbsorptionModel"));
    assert(document.at("Pos").contains("SrcDepth"));
    assert(document.at("Pos").contains("RecvRange"));
    assert(document.at("Pos").contains("RecvDepth"));
    assert(document.at("Pos").contains("RecvAzim"));
    assert(!document.at("Pos").contains("Sz"));
    assert(document.at("freqinfo").at("freq").is_number());
    assert(document.at("sspInput").front().at("HSTop").at("BC") == "vacuum");
    assert(document.at("ReflectionCoef").at("RTop").contains("R"));
    assert(!document.at("ReflectionCoef").contains("isDeg"));
    assert(document.at("SBP").is_object() || document.at("SBP").is_null());
    assert(source.to_json(output.string()));

    Interface restored;
    assert(restored.from_json(output.string()));
    const auto &a = source.getParams_const();
    const auto &b = restored.getParams_const();
    assert(a.Title == b.Title);
    assert(a.freqinfo.freq == b.freqinfo.freq);
    assert(a.Pos.NSz == b.Pos.NSz);
    assert(a.Pos.NRr == b.Pos.NRr);
    assert(a.sspInput.size() == b.sspInput.size());
    assert(a.sspInput.front().layers.front().alphaR.isApprox(
        b.sspInput.front().layers.front().alphaR));
    assert(a.sspInput.front().enableRootRestarts ==
           b.sspInput.front().enableRootRestarts);
    assert(a.SBP.NSBPPts == b.SBP.NSBPPts);
    assert(b.is_Velocity);

    Interface setters;
    setters.set_Title("setter json");
    setters.set_Freq(75.0);
    setters.set_Sz(10.0, 30.0, 3);
    setters.set_Rz((Eigen::VectorXd(2) << 20.0, 40.0).finished());
    setters.set_Rr(1000.0, 3000.0, 3);
    setters.set_Ro(Eigen::VectorXd::Zero(2));
    const OpenOcean_json setterDocument = OpenOcean_json::parse(setters.to_json_string());
    OOKC_parameters setterRoundTrip;
    from_json(setterDocument, setterRoundTrip);
    assert(setterRoundTrip.Title == "setter json");
    assert(setterRoundTrip.Pos.NSz == 3);
    assert(setterRoundTrip.Pos.is_Linspace_Sz);
    assert(setterRoundTrip.Pos.NRr == 3);

    OpenOcean_json aliasDocument = document;
    aliasDocument.erase("schemaVersion");
    OOKC_parameters aliasRoundTrip;
    from_json(aliasDocument, aliasRoundTrip);
    assert(aliasRoundTrip.Title == source.getParams_const().Title);

    OpenOcean_json multiFrequency = document;
    multiFrequency["freqinfo"]["Nfreq"] = 2;
    multiFrequency["freqinfo"]["freq"] = OpenOcean_json::array({50.0, 100.0});
    OOKC_parameters multiFrequencyRoundTrip;
    from_json(multiFrequency, multiFrequencyRoundTrip);
    assert(multiFrequencyRoundTrip.freqinfo.Nfreq == 2);
    assert(multiFrequencyRoundTrip.freqinfo.freqvec.size() == 2);
    assert(multiFrequencyRoundTrip.freqinfo.freqvec[1] == 100.0);

    OpenOcean_json capabilityRoundTripDocument = document;
    capabilityRoundTripDocument["AttenUnit"]["OceanAbsorptionModel"] = "Thorpe";
    capabilityRoundTripDocument["sspInput"].front()["SSPType"] = "cPCHIP";
    capabilityRoundTripDocument["sspInput"].front()["HSTop"]["BC"] = "grain";
    OOKC_parameters capabilityRoundTrip;
    from_json(capabilityRoundTripDocument, capabilityRoundTrip);
    assert(capabilityRoundTrip.AttenUnit.absModel == OceanAbsorptionModel::Thorpe);
    assert(capabilityRoundTrip.sspInput.front().SSPType == SSP_Mode::MODE_P_cPCHIP);
    assert(capabilityRoundTrip.sspInput.front().HSTop.BC == BC_Mode::MODE_G_Grain);
    OpenOcean_json capabilityReencoded;
    to_json(capabilityReencoded, capabilityRoundTrip);
    assert(capabilityReencoded["AttenUnit"]["OceanAbsorptionModel"] == "Thorpe");
    assert(capabilityReencoded["sspInput"].front()["SSPType"] == "cPCHIP");
    assert(capabilityReencoded["sspInput"].front()["HSTop"]["BC"] == "grain");

    OpenOcean_json invalidMultiFrequency = multiFrequency;
    invalidMultiFrequency["sspInput"].front()["HSBot"]["rho"] = -1.0;
    bool rejectedInvalidMultiFrequency = false;
    try
    {
        OOKC_parameters ignoredParams;
        from_json(invalidMultiFrequency, ignoredParams);
    }
    catch (const std::exception &)
    {
        rejectedInvalidMultiFrequency = true;
    }
    assert(rejectedInvalidMultiFrequency);

    {
        std::ofstream stream(invalid);
        stream << R"({"schemaVersion":1,"Title":"bad"})";
    }
    const std::string oldTitle = restored.getParams_const().Title;
    assert(!restored.from_json(invalid.string()));
    assert(restored.getParams_const().Title == oldTitle);

    OpenOcean_json invalidPhysical = document;
    invalidPhysical["Pos"]["RecvRange"] = OpenOcean_json::array();
    invalidPhysical.erase("paths");
    {
        std::ofstream stream(invalid);
        stream << invalidPhysical.dump(4);
    }
    assert(!restored.from_json(invalid.string()));
    assert(restored.getParams_const().Title == oldTitle);

    std::error_code ignored;
    std::filesystem::remove(output, ignored);
    std::filesystem::remove(invalid, ignored);
    return 0;
}
