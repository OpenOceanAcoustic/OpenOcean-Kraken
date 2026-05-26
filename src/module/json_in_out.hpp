#ifndef Kraken_JSON_IN_OUT_HPP
#define Kraken_JSON_IN_OUT_HPP

#include "json_eigen.hpp"
#include "OpenOceanKrakenParams.h"
#include "OpenOceanKrakenInterface.h"

// 定义数据结构对象的值
namespace OpenOceanKraken
{

    // 各类enum
    void to_json(OpenOcean_json &out, const SSP_Mode &mode)
    {
        switch (mode)
        {
        case SSP_Mode::MODE_A_Analytic:
            out = "Analytic";
            break;
        case SSP_Mode::MODE_C_cLinear:
            out = "cLinear";
            break;
        case SSP_Mode::MODE_N_n2Linear:
            out = "n2Linear";
            break;
        case SSP_Mode::MODE_P_cPCHIP:
            out = "cPCHIP";
            break;
        case SSP_Mode::MODE_S_cCubic:
            out = "cCubic";
            break;
        default:
            throw std::runtime_error("Unknown SSP_Mode value");
        }
    }
    
    void from_json(const OpenOcean_json &in, SSP_Mode &mode)
    {
        if (!in.is_string())
        {
            throw std::runtime_error("SSP_Mode must be a JSON string (e.g., \"Analytic\", \"cLinear\", etc.)");
        }

        const std::string &str = in.get_ref<const std::string &>();

        if (str == "Analytic")
        {
            mode = SSP_Mode::MODE_A_Analytic;
        }
        else if (str == "cLinear")
        {
            mode = SSP_Mode::MODE_C_cLinear;
        }
        else if (str == "n2Linear")
        {
            mode = SSP_Mode::MODE_N_n2Linear;
        }
        else if (str == "cPCHIP")
        {
            mode = SSP_Mode::MODE_P_cPCHIP;
        }
        else if (str == "cCubic")
        {
            mode = SSP_Mode::MODE_S_cCubic;
        }
        else
        {
            throw std::runtime_error("Invalid SSP_Mode string: \"" + str + "\"");
        }
    }

    void to_json(OpenOcean_json &out, const Media_Mode &mode)
    {
        switch (mode)
        {
        case Media_Mode::MODE_A_Acoustic:
            out = "Acoustic";
            break;
        case Media_Mode::MODE_E_Elastic:
            out = "Elastic";
            break;
        default:
            throw std::runtime_error("Unknown Media_Mode value");
        }
    }
    void from_json(const OpenOcean_json &in, Media_Mode &mode)
    {
        if (!in.is_string())
        {
            throw std::runtime_error("Media_Mode must be a JSON string (e.g., \"Acoustic\", \"Elastic\", etc.)");
        }

        const std::string &str = in.get_ref<const std::string &>();

        if (str == "Acoustic")
        {
            mode = Media_Mode::MODE_A_Acoustic;
        }
        else if (str == "Elastic")
        {
            mode = Media_Mode::MODE_E_Elastic;
        }
        else
        {
            throw std::runtime_error("Invalid Media_Mode string: \"" + str + "\"");
        }
    }

    void to_json(OpenOcean_json &out, const AttenuationUnit &mode)
    {
        switch (mode)
        {
        case AttenuationUnit::MODE_F_dB_per_m_kHz:
            out = "dB/m/kHz";
            break;
        case AttenuationUnit::MODE_L_params_lose:
            out = "Params_Lose";
            break;
        case AttenuationUnit::MODE_M_dB_per_m:
            out = "dB/m";
            break;
        case AttenuationUnit::MODE_N_Nepers_per_m:
            out = "Nepers/m";
            break;
        case AttenuationUnit::MODE_Q_Quality_Factor:
            out = "Quality_Factor";
            break;
        case AttenuationUnit::MODE_W_db_per_lambda:
            out = "dB/lambda";
            break;

        default:
            throw std::runtime_error("Unknown AttenuationUnit value");
        }
    }

    void from_json(const OpenOcean_json &in, AttenuationUnit &mode)
    {
        if (!in.is_string())
        {
            throw std::runtime_error("AttenuationUnit must be a JSON string (e.g., \"dB/m/kHz\", \"Params_Lose\", etc.)");
        }
        const std::string &str = in.get_ref<const std::string &>();
        if (str == "dB/m/kHz")
        {
            mode = AttenuationUnit::MODE_F_dB_per_m_kHz;
        }
        else if (str == "Params_Lose")
        {
            mode = AttenuationUnit::MODE_L_params_lose;
        }
        else if (str == "dB/m")
        {
            mode = AttenuationUnit::MODE_M_dB_per_m;
        }
        else if (str == "Nepers/m")
        {
            mode = AttenuationUnit::MODE_N_Nepers_per_m;
        }
        else if (str == "Quality_Factor")
        {
            mode = AttenuationUnit::MODE_Q_Quality_Factor;
        }
        else if (str == "dB/lambda")
        {
            mode = AttenuationUnit::MODE_W_db_per_lambda;
        }
        else
        {
            throw std::runtime_error("Invalid AttenuationUnit string: \"" + str + "\"");
        }
    }

    void to_json(OpenOcean_json &out, const OceanAbsorptionModel &mode)
    {
        switch (mode)
        {
        case OceanAbsorptionModel::FrancGarr:
            out = "FrancGarr";
            break;
        case OceanAbsorptionModel::Thorpe:
            out = "Thorpe";
            break;
        case OceanAbsorptionModel::None:
            out = "None";
            break;
        default:
            throw std::runtime_error("Unknown OceanAbsorptionModel value");
        }
    }

    void from_json(const OpenOcean_json &in, OceanAbsorptionModel &mode)
    {
        if (!in.is_string())
        {
            throw std::runtime_error("OceanAbsorptionModel must be a JSON string (e.g., \"FrancGarr\", \"Thorpe\", \"None\", etc.)");
        }
        const std::string &str = in.get_ref<const std::string &>();
        if (str == "FrancGarr")
        {
            mode = OceanAbsorptionModel::FrancGarr;
        }
        else if (str == "Thorpe")
        {
            mode = OceanAbsorptionModel::Thorpe;
        }
        else if (str == "None")
        {
            mode = OceanAbsorptionModel::None;
        }
        else
        {
            throw std::runtime_error("Invalid OceanAbsorptionModel string: \"" + str + "\"");
        }
    }

    void to_json(OpenOcean_json &out, const Atten_Mode &mode)
    {
        out = OpenOcean_json{
            {"AttenuationUnit", mode.attnUnit},
            {"OceanAbsorptionModel", mode.absModel}};
    }
    void from_json(const OpenOcean_json &in, Atten_Mode &mode)
    {
        if (!in.is_object())
        {
            throw std::runtime_error("Atten_Mode must be a JSON object with fields \"AttenuationUnit\" and \"OceanAbsorptionModel\".");
        }
        mode.attnUnit = in.at("AttenuationUnit").get<AttenuationUnit>();
        mode.absModel = in.at("OceanAbsorptionModel").get<OceanAbsorptionModel>();
    }

    void to_json(OpenOcean_json &out, const BC_Mode &mode)
    {
        switch (mode)
        {
        case BC_Mode::MODE_A_Half_space:
            out = "halfspace";
            break;
        case BC_Mode::MODE_F_File:
            out = "file";
            break;
        case BC_Mode::MODE_G_Grain:
            out = "grain";
            break;
        case BC_Mode::MODE_P_Precomputed:
            out = "precomputed";
            break;
        case BC_Mode::MODE_R_Rigid:
            out = "rigid";
            break;
        case BC_Mode::MODE_V_Vacuum:
            out = "vacuum";
            break;
        default:
            throw std::runtime_error("Unknown BC_Mode value");
        }
    }
    void from_json(const OpenOcean_json &in, BC_Mode &mode)
    {
        if (!in.is_string())
        {
            throw std::runtime_error("BC_Mode must be a JSON string (e.g., \"halfspace\", \"file\", \"grain\", \"precomputed\", \"rigid\", \"vacuum\", etc.)");
        }
        const std::string &str = in.get_ref<const std::string &>();
        if (str == "halfspace")
        {
            mode = BC_Mode::MODE_A_Half_space;
        }
        else if (str == "file")
        {
            mode = BC_Mode::MODE_F_File;
        }
        else if (str == "grain")
        {
            mode = BC_Mode::MODE_G_Grain;
        }
        else if (str == "precomputed")
        {
            mode = BC_Mode::MODE_P_Precomputed;
        }
        else if (str == "rigid")
        {
            mode = BC_Mode::MODE_R_Rigid;
        }
        else if (str == "vacuum")
        {
            mode = BC_Mode::MODE_V_Vacuum;
        }
        else
        {
            throw std::runtime_error("Invalid BC_Mode string: \"" + str + "\"");
        }
    }

    void to_json(OpenOcean_json &out, const Source_Mode &mode)
    {
        switch (mode)
        {
        case Source_Mode::MODE_R_Point:
            out = "Point";
            break;
        case Source_Mode::MODE_X_Line:
            out = "Line";
            break;
        default:
            throw std::runtime_error("Unknown Source_Mode value");
        }
    }

    void from_json(const OpenOcean_json &in, Source_Mode &mode)
    {
        if (!in.is_string())
        {
            throw std::runtime_error("Source_Mode must be a JSON string (e.g., \"Point\", \"Line\", etc.)");
        }
        const std::string &str = in.get_ref<const std::string &>();
        if (str == "Point")
        {
            mode = Source_Mode::MODE_R_Point;
        }
        else if (str == "Line")
        {
            mode = Source_Mode::MODE_X_Line;
        }
        else
        {
            throw std::runtime_error("Invalid Source_Mode string: \"" + str + "\"");
        }
    }

    void to_json(OpenOcean_json &out, const Grid_Mode &mode)
    {
        switch (mode)
        {
        case Grid_Mode::MODE_I_Irregular:
            out = "Irregular";
            break;
        case Grid_Mode::MODE_R_Rectangular:
            out = "Rectangular";
            break;
        default:
            throw std::runtime_error("Unknown Grid_Mode value");
        }
    }
    void from_json(const OpenOcean_json &in, Grid_Mode &mode)
    {
        if (!in.is_string())
        {
            throw std::runtime_error("Grid_Mode must be a JSON string (e.g., \"Irregular\", \"Rectangular\", etc.)");
        }
        const std::string &str = in.get_ref<const std::string &>();
        if (str == "Irregular")
        {
            mode = Grid_Mode::MODE_I_Irregular;
        }
        else if (str == "Rectangular")
        {
            mode = Grid_Mode::MODE_R_Rectangular;
        }
        else
        {
            throw std::runtime_error("Invalid Grid_Mode string: \"" + str + "\"");
        }
    }

    void to_json(OpenOcean_json &out, const Run_Mode &mode)
    {
        switch (mode)
        {
        case Run_Mode::MODE_F_Field:
            out = "Field";
            break;
        case Run_Mode::MODE_M_Modes:
            out = "Modes";
            break;
        case Run_Mode::MODE_B_Both:
            out = "Both";
            break;

        default:
            throw std::runtime_error("Unknown Run_Mode value");
        }
    }
    void from_json(const OpenOcean_json &in, Run_Mode &mode)
    {
        if (!in.is_string())
        {
            throw std::runtime_error("Run_Mode must be a JSON string (e.g., \"Field\", \"Modes\", \"Both\", etc.)");
        }
        const std::string &str = in.get_ref<const std::string &>();
        if (str == "Field")
        {
            mode = Run_Mode::MODE_F_Field;
        }
        else if (str == "Modes")
        {
            mode = Run_Mode::MODE_M_Modes;
        }
        else if (str == "Both")
        {
            mode = Run_Mode::MODE_B_Both;
        }
        else
        {
            throw std::runtime_error("Invalid Run_Mode string: \"" + str + "\"");
        }
    }

    void to_json(OpenOcean_json &out, const CoherenceType &mode)
    {
        switch (mode)
        {
        case CoherenceType::Coherent:
            out = "Coherent";
            break;
        case CoherenceType::Incoherent:
            out = "Incoherent";
            break;
        default:
            throw std::runtime_error("Unknown CoherenceType value");
        }
    }
    void from_json(const OpenOcean_json &in, CoherenceType &mode)
    {
        if (!in.is_string())
        {
            throw std::runtime_error("CoherenceType must be a JSON string (e.g., \"Coherent\", \"Incoherent\", etc.)");
        }
        const std::string &str = in.get_ref<const std::string &>();
        if (str == "Coherent")
        {
            mode = CoherenceType::Coherent;
        }
        else if (str == "Incoherent")
        {
            mode = CoherenceType::Incoherent;
        }
        else
        {
            throw std::runtime_error("Invalid CoherenceType string: \"" + str + "\"");
        }
    }

    void to_json(OpenOcean_json &out, const ModeType &mode)
    {
        switch (mode)
        {
        case ModeType::Adiabatic:
            out = "Adiabatic";
            break;
        case ModeType::Couple:
            out = "Couple";
            break;
        default:
            throw std::runtime_error("Unknown ModeType value");
        }
    }

    void from_json(const OpenOcean_json &in, ModeType &mode)
    {
        if (!in.is_string())
        {
            throw std::runtime_error("ModeType must be a JSON string (e.g., \"Adiabatic\", \"Couple\", etc.)");
        }
        const std::string &str = in.get_ref<const std::string &>();
        if (str == "Adiabatic")
        {
            mode = ModeType::Adiabatic;
        }
        else if (str == "Couple")
        {
            mode = ModeType::Couple;
        }
        else
        {
            throw std::runtime_error("Invalid ModeType string: \"" + str + "\"");
        }
    }

    // HSInfo
    void to_json(OpenOcean_json &out, const HSInfo &HS)
    {
        out = OpenOcean_json{
            {"Depth", HS.Depth},
            {"alphaR", HS.alphaR},
            {"alphaI", HS.alphaI},
            {"betaR", HS.betaR},
            {"betaI", HS.betaI},
            {"rho", HS.rho},
            {"beta", HS.beta},
            {"ft", HS.ft},
            {"BC", HS.BC} // 假设 BC_Mode 已有 to_json
        };
    }
    void from_json(const OpenOcean_json &in, HSInfo &HS)
    {
        if (!in.is_object())
        {
            throw std::runtime_error("HSInfo must be a JSON object");
        }

        HS.Depth = in.at("Depth").get<double>();
        HS.alphaR = in.at("alphaR").get<double>();
        HS.alphaI = in.at("alphaI").get<double>();
        HS.betaR = in.at("betaR").get<double>();
        HS.betaI = in.at("betaI").get<double>();
        HS.rho = in.at("rho").get<double>();
        HS.beta = in.at("beta").get<double>();
        HS.ft = in.at("ft").get<double>();
        HS.BC = in.at("BC").get<BC_Mode>(); // 自动调用 BC_Mode 的 from_json
    }

    // position
    void to_json(OpenOcean_json &out, const Position &pos)
    {
        OpenOcean_json j;

        j["GridType"] = pos.GridType;

        // === Sz ===
        if (pos.is_Linspace_Sz && pos.Sz.size() >= 2)
        {
            j["SrcDepth"]["start"] = pos.Sz(0);
            j["SrcDepth"]["end"] = pos.Sz(pos.Sz.size() - 1);
            j["SrcDepth"]["NSz"] = pos.NSz;
        }
        else
        {
            j["SrcDepth"] = pos.Sz;
        }

        // === Rr ===
        if (pos.is_Linspace_Rr && pos.Rr.size() >= 2)
        {
            j["RecvRange"]["start"] = pos.Rr(0);
            j["RecvRange"]["end"] = pos.Rr(pos.Rr.size() - 1);
            j["RecvRange"]["NRr"] = pos.NRr;
        }
        else
        {
            j["RecvRange"] = pos.Rr;
        }

        // === Rz ===
        if (pos.is_Linspace_Rz && pos.Rz.size() >= 2)
        {
            j["RecvDepth"]["start"] = pos.Rz(0);
            j["RecvDepth"]["end"] = pos.Rz(pos.Rz.size() - 1);
            j["RecvDepth"]["NRz"] = pos.NRz;
        }
        else
        {
            j["RecvDepth"] = pos.Rz;
        }

        // === Ro ===
        if (pos.is_Linspace_Ro && pos.Ro.size() >= 2)
        {
            j["RecvAzim"]["start"] = pos.Ro(0);
            j["RecvAzim"]["end"] = pos.Ro(pos.Ro.size() - 1);
            j["RecvAzim"]["NRo"] = pos.NRo;
        }
        else
        {
            j["RecvAzim"] = pos.Ro;
        }

        out = std::move(j);
    }
    void from_json(const OpenOcean_json &in, Position &pos)
    {
        // GridType
        pos.GridType = in.at("GridType").get<Grid_Mode>();

        // === SrcDepth → Sz ===
        const auto &src_depth = in.at("SrcDepth");
        if (src_depth.is_object() && src_depth.contains("start"))
        {
            double start = src_depth.at("start").get<double>();
            double end = src_depth.at("end").get<double>();
            int N = src_depth.at("NSz").get<int>();
            pos.Sz = Eigen::VectorXd::LinSpaced(N, start, end);
            pos.is_Linspace_Sz = true;
        }
        else
        {
            pos.Sz = src_depth.get<Eigen::VectorXd>();
            pos.is_Linspace_Sz = false;
        }

        // === RecvRange → Rr ===
        const auto &recv_range = in.at("RecvRange");
        if (recv_range.is_object() && recv_range.contains("start"))
        {
            double start = recv_range.at("start").get<double>();
            double end = recv_range.at("end").get<double>();
            int N = recv_range.at("NRr").get<int>();
            pos.Rr = Eigen::VectorXd::LinSpaced(N, start, end);
            pos.is_Linspace_Rr = true;
        }
        else
        {
            pos.Rr = recv_range.get<Eigen::VectorXd>();
            pos.is_Linspace_Rr = false;
        }

        // === RecvDepth → Rz ===
        const auto &recv_depth = in.at("RecvDepth");
        if (recv_depth.is_object() && recv_depth.contains("start"))
        {
            double start = recv_depth.at("start").get<double>();
            double end = recv_depth.at("end").get<double>();
            int N = recv_depth.at("NRz").get<int>();
            pos.Rz = Eigen::VectorXd::LinSpaced(N, start, end);
            pos.is_Linspace_Rz = true;
        }
        else
        {
            pos.Rz = recv_depth.get<Eigen::VectorXd>();
            pos.is_Linspace_Rz = false;
        }

        // === RecvAzim → Ro ===
        const auto &recv_azim = in.at("RecvAzim");
        if (recv_azim.is_object() && recv_azim.contains("start"))
        {
            double start = recv_azim.at("start").get<double>();
            double end = recv_azim.at("end").get<double>();
            int N = recv_azim.at("NRo").get<int>();
            pos.Ro = Eigen::VectorXd::LinSpaced(N, start, end);
            pos.is_Linspace_Ro = true;
        }
        else
        {
            pos.Ro = recv_azim.get<Eigen::VectorXd>();
            pos.is_Linspace_Ro = false;
        }
    }

    // Reflection
    // void to_json(OpenOcean_json &out, const ReflectionCoef &ref)
    // {
    //     out = {{"theta", ref.theta}, {"R", ref.R}, {"phi", ref.phi}};
    // }
    // void from_json(const OpenOcean_json &in, ReflectionCoef &ref)
    // {
    //     ref.theta = in.at("theta").get<double>();
    //     ref.R = in.at("R").get<double>();
    //     ref.phi = in.at("phi").get<double>();
    // }
    // ReflectionCoefInfo
    void to_json(OpenOcean_json &out, const ReflectionCoefInfo &ref)
    {
        int TOP_Size = ref.RTop.size();
        int BOT_Size = ref.RBot.size();

        // 转换
        Eigen::VectorXd RTop_R(TOP_Size);
        Eigen::VectorXd RBot_R(BOT_Size);
        Eigen::VectorXd RTop_phi(TOP_Size);
        Eigen::VectorXd RBot_phi(BOT_Size);
        Eigen::VectorXd RTot_theta(TOP_Size);
        Eigen::VectorXd RBot_theta(BOT_Size);
        double factor = ref.isDeg ? RadDeg : 1.0;
        for (int i = 0; i < TOP_Size; i++)
        {
            RTop_R(i) = ref.RTop[i].R;
            RTop_phi(i) = ref.RTop[i].phi * factor;
            RTot_theta(i) = ref.RTop[i].theta;
        }
        for (int i = 0; i < BOT_Size; i++)
        {
            RBot_R(i) = ref.RBot[i].R;
            RBot_phi(i) = ref.RBot[i].phi * factor;
            RBot_theta(i) = ref.RBot[i].theta;
        }
        OpenOcean_json j;
        j["RTop"]["R"] = RTop_R;
        j["RTop"]["phi"] = RTop_phi;
        j["RTop"]["theta"] = RTot_theta;
        j["RBot"]["R"] = RBot_R;
        j["RBot"]["phi"] = RBot_phi;
        j["RBot"]["theta"] = RBot_theta;
        out = std::move(j);
    }
    void from_json(const OpenOcean_json &in, ReflectionCoefInfo &ref)
    {
        auto load_coefs = [&](const OpenOcean_json &obj) -> Eigen::Matrix<ReflectionCoef, 1, Eigen::Dynamic>
        {
            const auto &R_vec = obj.at("R").get<Eigen::VectorXd>();
            const auto &phi_vec = obj.at("phi").get<Eigen::VectorXd>();
            const auto &theta_vec = obj.at("theta").get<Eigen::VectorXd>();

            int N = static_cast<int>(R_vec.size());
            Eigen::Matrix<ReflectionCoef, 1, Eigen::Dynamic> coefs(1, N);

            for (int i = 0; i < N; ++i)
            {
                coefs(i).R = R_vec(i);
                coefs(i).phi = phi_vec(i);
                coefs(i).theta = theta_vec(i);
            }
            return coefs;
        };

        ref.RTop = load_coefs(in.at("RTop"));
        ref.RBot = load_coefs(in.at("RBot"));
        ref.isDeg = false;
    }

    // SrcBmPat
    void to_json(OpenOcean_json &out, const SrcBmPat &pat)
    {
        if (pat.isSet)
        {
            out["NSBPPts"] = pat.NSBPPts;
            out["theta"] = pat.theta;
            out["pat"] = pat.pat;
        }
    }
    void from_json(const OpenOcean_json &in, SrcBmPat &pat)
    {

            pat.NSBPPts = in.at("NSBPPts").get<int>();
            pat.theta = in.at("theta").get<Eigen::VectorXd>();
            pat.pat = in.at("pat").get<Eigen::VectorXd>();
        
        pat.isSet = true;
    }
    // freqInfo
    void to_json(OpenOcean_json &out, const FreqInfo &freq)
    {
        if (freq.Nfreq > 1)
        {
            out["Nfreq"] = freq.Nfreq;
            out["freq"] = freq.freqvec;
        }
        else
        {
            out["Nfreq"] = freq.Nfreq;
            out["freq"] = freq.freq;
        }
    }
    void from_json(const OpenOcean_json &in, FreqInfo &freq)
    {

            freq.Nfreq = in.at("Nfreq").get<int>();
            if (freq.Nfreq > 1)
            {
                freq.freqvec = in.at("freq").get<Eigen::VectorXd>();
            }
            else
            {
                freq.freq = in.at("freq").get<double>();
            }
        
    }

    // ssp命名空间里面的转换
    namespace ssp
    {

        // SSPLayer
        void to_json(OpenOcean_json &out, const SSPLayer &Layer)
        {
            OpenOcean_json j;
            j["npts"] = Layer.npts;
            j["nmesh"] = Layer.nmesh;
            j["beta"] = Layer.beta;
            j["ft"] = Layer.ft;
            j["sigma"] = Layer.sigma;
            j["z"] = Layer.z;
            j["rho"] = Layer.rho;
            j["alphaR"] = Layer.alphaR;
            j["alphaI"] = Layer.alphaI;
            j["betaR"] = Layer.betaR;
            j["betaI"] = Layer.betaI;
            j["Material"] = Layer.Material;
            out = std::move(j);
        }
        void from_json(const OpenOcean_json &in, SSPLayer &Layer)
        {

            Layer.npts = in.at("npts").get<int>();
            Layer.nmesh = in.at("nmesh").get<int>();
            Layer.beta = in.at("beta").get<double>();
            Layer.ft = in.at("ft").get<double>();
            Layer.sigma = in.at("sigma").get<double>();
            Layer.z = in.at("z").get<Eigen::VectorXd>();
            Layer.rho = in.at("rho").get<Eigen::VectorXd>();
            Layer.alphaR = in.at("alphaR").get<Eigen::VectorXd>();
            Layer.alphaI = in.at("alphaI").get<Eigen::VectorXd>();
            Layer.betaR = in.at("betaR").get<Eigen::VectorXd>();
            Layer.betaI = in.at("betaI").get<Eigen::VectorXd>();
            Layer.Material = in.at("Material").get<Media_Mode>();
        }

        // Range_Independent_Area
        void to_json(OpenOcean_json &out, const Range_Independent_Area &Area)
        {
            OpenOcean_json j;
            j["SSPType"] = Area.SSPType;
            j["Range"] = Area.Range;
            j["layers"] = Area.layers;
            j["HSTop"] = Area.HSTop;
            j["HSBot"] = Area.HSBot;
            out = std::move(j);
        }
        void from_json(const OpenOcean_json &in, Range_Independent_Area &Area)
        {
            Area.SSPType = in.at("SSPType").get<SSP_Mode>();
            Area.Range = in.at("Range").get<double>();
            Area.layers = in.at("layers").get<std::vector<SSPLayer>>();
            Area.HSTop = in.at("HSTop").get<HSInfo>();
            Area.HSBot = in.at("HSBot").get<HSInfo>();
        }
    }

    // params
    void to_json(OpenOcean_json &out, const OOK_parameters &params)
    {
        OpenOcean_json j;
        j["Title"] = params.Title;
        j["freqinfo"] = params.freqinfo;
        j["AttenUnit"] = params.AttenUnit;
        j["Pos"] = params.Pos;
        j["sspInput"] = params.sspInput;
        j["ReflectionCoef"] = params.ReflectionCoef;
        j["SBP"] = params.SBP;
        j["is_Velocity"] = params.is_Velocity;
        j["cLow"] = params.cLow;
        j["cHigh"] = params.cHigh;
        j["Rmax"] = params.Rmax;
        j["SourceType"] = params.SourceType;
        j["RunMode"] = params.runMode;
        j["CoherenceType"] = params.coherenceType;
        j["ModeType"] = params.modeType;
        out = std::move(j);
    }

    void from_json(const OpenOcean_json &in, OOK_parameters &params)
    {
        params.Title = in.at("Title").get<std::string>();
        params.freqinfo = in.at("freqinfo").get<FreqInfo>();
        params.AttenUnit = in.at("AttenUnit").get<Atten_Mode>();
        params.Pos = in.at("Pos").get<Position>();
        params.sspInput = in.at("sspInput").get<std::vector<ssp::Range_Independent_Area>>();
        params.ReflectionCoef = in.at("ReflectionCoef").get<ReflectionCoefInfo>();
        params.SBP = in.at("SBP").get<SrcBmPat>();
        params.is_Velocity = in.at("is_Velocity").get<bool>();
        params.cLow = in.at("cLow").get<double>();
        params.cHigh = in.at("cHigh").get<double>();
        params.Rmax = in.at("Rmax").get<double>();
        params.SourceType = in.at("SourceType").get<Source_Mode>();
        params.runMode = in.at("RunMode").get<Run_Mode>();
        params.coherenceType = in.at("CoherenceType").get<CoherenceType>();
        params.modeType = in.at("ModeType").get<ModeType>();
        
    }

}

#endif
