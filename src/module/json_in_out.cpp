#include "module/json_in_out.hpp"
#include "module/ParameterAdapters.h"

#include <cmath>
#include <complex>
#include <fstream>
#include <stdexcept>
#include <utility>

namespace OpenOceanKrakenc
{
namespace
{
template <typename Enum>
void enumToJson(OpenOcean_json &out, Enum value,
                std::initializer_list<std::pair<Enum, const char *>> names,
                const char *typeName)
{
    for (const auto &entry : names)
    {
        if (entry.first == value)
        {
            out = entry.second;
            return;
        }
    }
    throw std::runtime_error(std::string("unknown ") + typeName + " value");
}

template <typename Enum>
void enumFromJson(const OpenOcean_json &in, Enum &value,
                  std::initializer_list<std::pair<Enum, const char *>> names,
                  const char *typeName)
{
    if (!in.is_string())
    {
        throw std::runtime_error(std::string(typeName) + " must be a JSON string");
    }
    const std::string text = in.get<std::string>();
    for (const auto &entry : names)
    {
        if (text == entry.second)
        {
            value = entry.first;
            return;
        }
    }
    throw std::runtime_error(std::string("invalid ") + typeName + ": " + text);
}

OpenOcean_json realVector(const Eigen::VectorXd &values)
{
    OpenOcean_json out = OpenOcean_json::array();
    for (Eigen::Index i = 0; i < values.size(); ++i) out.push_back(values[i]);
    return out;
}

OpenOcean_json intVector(const Eigen::VectorXi &values)
{
    OpenOcean_json out = OpenOcean_json::array();
    for (Eigen::Index i = 0; i < values.size(); ++i) out.push_back(values[i]);
    return out;
}

OpenOcean_json complexValue(const std::complex<double> &value)
{
    return OpenOcean_json::array({value.real(), value.imag()});
}

OpenOcean_json complexVector(const Eigen::VectorXcd &values)
{
    OpenOcean_json out = OpenOcean_json::array();
    for (Eigen::Index i = 0; i < values.size(); ++i) out.push_back(complexValue(values[i]));
    return out;
}

OpenOcean_json complexMatrix(const Eigen::MatrixXcd &values)
{
    OpenOcean_json data = OpenOcean_json::array();
    for (Eigen::Index row = 0; row < values.rows(); ++row)
    {
        OpenOcean_json line = OpenOcean_json::array();
        for (Eigen::Index col = 0; col < values.cols(); ++col)
        {
            line.push_back(complexValue(values(row, col)));
        }
        data.push_back(std::move(line));
    }
    return data;
}

Eigen::VectorXd readRealVector(const OpenOcean_json &in, const char *name)
{
    if (!in.is_array()) throw std::runtime_error(std::string(name) + " must be an array");
    Eigen::VectorXd out(static_cast<Eigen::Index>(in.size()));
    for (Eigen::Index i = 0; i < out.size(); ++i) out[i] = in.at(static_cast<std::size_t>(i)).get<double>();
    return out;
}

Eigen::VectorXi readIntVector(const OpenOcean_json &in, const char *name)
{
    if (!in.is_array()) throw std::runtime_error(std::string(name) + " must be an array");
    Eigen::VectorXi out(static_cast<Eigen::Index>(in.size()));
    for (Eigen::Index i = 0; i < out.size(); ++i) out[i] = in.at(static_cast<std::size_t>(i)).get<int>();
    return out;
}

std::complex<double> readComplex(const OpenOcean_json &in, const char *name)
{
    if (!in.is_array() || in.size() != 2)
    {
        throw std::runtime_error(std::string(name) + " must be [real, imaginary]");
    }
    return {in.at(0).get<double>(), in.at(1).get<double>()};
}

Eigen::VectorXcd readComplexVector(const OpenOcean_json &in, const char *name)
{
    if (!in.is_array()) throw std::runtime_error(std::string(name) + " must be an array");
    Eigen::VectorXcd out(static_cast<Eigen::Index>(in.size()));
    for (Eigen::Index i = 0; i < out.size(); ++i) out[i] = readComplex(in.at(static_cast<std::size_t>(i)), name);
    return out;
}

Eigen::MatrixXcd readComplexMatrix(const OpenOcean_json &in, const char *name)
{
    if (!in.is_array()) throw std::runtime_error(std::string(name) + " must be a matrix");
    const Eigen::Index rows = static_cast<Eigen::Index>(in.size());
    const Eigen::Index cols = rows == 0 ? 0 : static_cast<Eigen::Index>(in.at(0).size());
    Eigen::MatrixXcd out(rows, cols);
    for (Eigen::Index row = 0; row < rows; ++row)
    {
        if (!in.at(static_cast<std::size_t>(row)).is_array() ||
            static_cast<Eigen::Index>(in.at(static_cast<std::size_t>(row)).size()) != cols)
        {
            throw std::runtime_error(std::string(name) + " has ragged rows");
        }
        for (Eigen::Index col = 0; col < cols; ++col)
        {
            out(row, col) = readComplex(in.at(static_cast<std::size_t>(row)).at(static_cast<std::size_t>(col)), name);
        }
    }
    return out;
}

void checkPosition(const Position &value)
{
    if (value.NSz != value.Sz.size() || value.NRr != value.Rr.size() ||
        value.NRz != value.Rz.size() || value.NRo != value.Ro.size())
    {
        throw std::runtime_error("Position counts must match vector sizes");
    }
    if (!value.Sz.allFinite() || !value.Rr.allFinite() ||
        !value.Rz.allFinite() || !value.Ro.allFinite())
    {
        throw std::runtime_error("Position vectors must contain finite values");
    }
    if (value.Ro.size() != 0 && value.Rz.size() != 0 &&
        value.Ro.size() != 1 && value.Ro.size() != value.Rz.size())
    {
        throw std::runtime_error(
            "Position RecvAzim must contain one value or match RecvDepth");
    }
}
}

void to_json(OpenOcean_json &out, const SSP_Mode &v) { enumToJson(out, v, {{SSP_Mode::MODE_N_n2Linear,"n2Linear"},{SSP_Mode::MODE_C_cLinear,"cLinear"},{SSP_Mode::MODE_P_cPCHIP,"cPCHIP"},{SSP_Mode::MODE_S_cCubic,"cCubic"},{SSP_Mode::MODE_A_Analytic,"Analytic"}}, "SSP_Mode"); }
void from_json(const OpenOcean_json &in, SSP_Mode &v) { enumFromJson(in, v, {{SSP_Mode::MODE_N_n2Linear,"n2Linear"},{SSP_Mode::MODE_C_cLinear,"cLinear"},{SSP_Mode::MODE_P_cPCHIP,"cPCHIP"},{SSP_Mode::MODE_S_cCubic,"cCubic"},{SSP_Mode::MODE_A_Analytic,"Analytic"}}, "SSP_Mode"); }
void to_json(OpenOcean_json &out, const Media_Mode &v) { enumToJson(out, v, {{Media_Mode::MODE_A_Acoustic,"Acoustic"},{Media_Mode::MODE_E_Elastic,"Elastic"}}, "Media_Mode"); }
void from_json(const OpenOcean_json &in, Media_Mode &v) { enumFromJson(in, v, {{Media_Mode::MODE_A_Acoustic,"Acoustic"},{Media_Mode::MODE_E_Elastic,"Elastic"}}, "Media_Mode"); }
void to_json(OpenOcean_json &out, const AttenuationUnit &v) { enumToJson(out, v, {{AttenuationUnit::MODE_F_dB_per_m_kHz,"dB/m/kHz"},{AttenuationUnit::MODE_L_params_lose,"Params_Lose"},{AttenuationUnit::MODE_M_dB_per_m,"dB/m"},{AttenuationUnit::MODE_m_dB_per_m,"dB/m-lowercase"},{AttenuationUnit::MODE_N_Nepers_per_m,"Nepers/m"},{AttenuationUnit::MODE_Q_Quality_Factor,"Quality_Factor"},{AttenuationUnit::MODE_W_db_per_lambda,"dB/lambda"}}, "AttenuationUnit"); }
void from_json(const OpenOcean_json &in, AttenuationUnit &v) { enumFromJson(in, v, {{AttenuationUnit::MODE_F_dB_per_m_kHz,"dB/m/kHz"},{AttenuationUnit::MODE_L_params_lose,"Params_Lose"},{AttenuationUnit::MODE_M_dB_per_m,"dB/m"},{AttenuationUnit::MODE_m_dB_per_m,"dB/m-lowercase"},{AttenuationUnit::MODE_N_Nepers_per_m,"Nepers/m"},{AttenuationUnit::MODE_Q_Quality_Factor,"Quality_Factor"},{AttenuationUnit::MODE_W_db_per_lambda,"dB/lambda"}}, "AttenuationUnit"); }
void to_json(OpenOcean_json &out, const OceanAbsorptionModel &v) { enumToJson(out, v, {{OceanAbsorptionModel::None,"None"},{OceanAbsorptionModel::Thorpe,"Thorpe"},{OceanAbsorptionModel::FrancGarr,"FrancGarr"}}, "OceanAbsorptionModel"); }
void from_json(const OpenOcean_json &in, OceanAbsorptionModel &v) { enumFromJson(in, v, {{OceanAbsorptionModel::None,"None"},{OceanAbsorptionModel::Thorpe,"Thorpe"},{OceanAbsorptionModel::FrancGarr,"FrancGarr"}}, "OceanAbsorptionModel"); }
void to_json(OpenOcean_json &out, const BC_Mode &v) { enumToJson(out, v, {{BC_Mode::MODE_R_Rigid,"rigid"},{BC_Mode::MODE_V_Vacuum,"vacuum"},{BC_Mode::MODE_F_File,"file"},{BC_Mode::MODE_A_Half_space,"halfspace"},{BC_Mode::MODE_G_Grain,"grain"},{BC_Mode::MODE_P_Precomputed,"precomputed"}}, "BC_Mode"); }
void from_json(const OpenOcean_json &in, BC_Mode &v) { enumFromJson(in, v, {{BC_Mode::MODE_R_Rigid,"rigid"},{BC_Mode::MODE_V_Vacuum,"vacuum"},{BC_Mode::MODE_F_File,"file"},{BC_Mode::MODE_A_Half_space,"halfspace"},{BC_Mode::MODE_G_Grain,"grain"},{BC_Mode::MODE_P_Precomputed,"precomputed"},{BC_Mode::MODE_R_Rigid,"Rigid"},{BC_Mode::MODE_V_Vacuum,"Vacuum"},{BC_Mode::MODE_F_File,"File"},{BC_Mode::MODE_A_Half_space,"HalfSpace"},{BC_Mode::MODE_G_Grain,"Grain"},{BC_Mode::MODE_P_Precomputed,"Precomputed"}}, "BC_Mode"); }
void to_json(OpenOcean_json &out, const Source_Mode &v) { enumToJson(out, v, {{Source_Mode::MODE_R_Point,"Point"},{Source_Mode::MODE_X_Line,"Line"},{Source_Mode::MODE_S_ScaledCylindrical,"ScaledCylindrical"}}, "Source_Mode"); }
void from_json(const OpenOcean_json &in, Source_Mode &v) { enumFromJson(in, v, {{Source_Mode::MODE_R_Point,"Point"},{Source_Mode::MODE_X_Line,"Line"},{Source_Mode::MODE_S_ScaledCylindrical,"ScaledCylindrical"}}, "Source_Mode"); }
void to_json(OpenOcean_json &out, const Grid_Mode &v) { enumToJson(out, v, {{Grid_Mode::MODE_R_Rectangular,"Rectangular"},{Grid_Mode::MODE_I_Irregular,"Irregular"}}, "Grid_Mode"); }
void from_json(const OpenOcean_json &in, Grid_Mode &v) { enumFromJson(in, v, {{Grid_Mode::MODE_R_Rectangular,"Rectangular"},{Grid_Mode::MODE_I_Irregular,"Irregular"}}, "Grid_Mode"); }
void to_json(OpenOcean_json &out, const Run_Mode &v) { enumToJson(out, v, {{Run_Mode::MODE_M_Modes,"Modes"},{Run_Mode::MODE_F_Field,"Field"},{Run_Mode::MODE_B_Both,"Both"}}, "Run_Mode"); }
void from_json(const OpenOcean_json &in, Run_Mode &v) { enumFromJson(in, v, {{Run_Mode::MODE_M_Modes,"Modes"},{Run_Mode::MODE_F_Field,"Field"},{Run_Mode::MODE_B_Both,"Both"}}, "Run_Mode"); }
void to_json(OpenOcean_json &out, const CoherenceType &v) { enumToJson(out, v, {{CoherenceType::Coherent,"Coherent"},{CoherenceType::Incoherent,"Incoherent"}}, "CoherenceType"); }
void from_json(const OpenOcean_json &in, CoherenceType &v) { enumFromJson(in, v, {{CoherenceType::Coherent,"Coherent"},{CoherenceType::Incoherent,"Incoherent"}}, "CoherenceType"); }
void to_json(OpenOcean_json &out, const ModeType &v) { enumToJson(out, v, {{ModeType::Adiabatic,"Adiabatic"},{ModeType::Couple,"Couple"}}, "ModeType"); }
void from_json(const OpenOcean_json &in, ModeType &v) { enumFromJson(in, v, {{ModeType::Adiabatic,"Adiabatic"},{ModeType::Couple,"Couple"}}, "ModeType"); }

void to_json(OpenOcean_json &out, const Atten_Mode &v) { out = {{"AttenuationUnit",v.attnUnit},{"OceanAbsorptionModel",v.absModel}}; }
void from_json(const OpenOcean_json &in, Atten_Mode &v) { const char *unit=in.contains("AttenuationUnit")?"AttenuationUnit":"attnUnit"; const char *model=in.contains("OceanAbsorptionModel")?"OceanAbsorptionModel":"absModel"; in.at(unit).get_to(v.attnUnit); in.at(model).get_to(v.absModel); }
void to_json(OpenOcean_json &out, const HSInfo &v) { out={{"alphaR",v.alphaR},{"alphaI",v.alphaI},{"betaR",v.betaR},{"betaI",v.betaI},{"beta",v.beta},{"ft",v.ft},{"cp",complexValue(v.cp)},{"cs",complexValue(v.cs)},{"rho",v.rho},{"Depth",v.Depth},{"BC",v.BC}}; }
void from_json(const OpenOcean_json &in, HSInfo &v) { v.alphaR=in.at("alphaR").get<double>(); v.alphaI=in.at("alphaI").get<double>(); v.betaR=in.at("betaR").get<double>(); v.betaI=in.at("betaI").get<double>(); v.beta=in.value("beta",0.0); v.ft=in.value("ft",0.0); v.cp=in.contains("cp")?readComplex(in.at("cp"),"cp"):std::complex<double>(v.alphaR,v.alphaI); v.cs=in.contains("cs")?readComplex(in.at("cs"),"cs"):std::complex<double>(v.betaR,v.betaI); v.rho=in.at("rho").get<double>(); v.Depth=in.at("Depth").get<double>(); in.at("BC").get_to(v.BC); }
void to_json(OpenOcean_json &out, const Position &v) { const auto encode=[](const Eigen::VectorXd &values,bool linspace,int count,const char *countName){ if(linspace&&values.size()>=2) return OpenOcean_json{{"start",values[0]},{"end",values[values.size()-1]},{countName,count}}; return realVector(values); }; out={{"GridType",v.GridType},{"SrcDepth",encode(v.Sz,v.is_Linspace_Sz,v.NSz,"NSz")},{"RecvRange",encode(v.Rr,v.is_Linspace_Rr,v.NRr,"NRr")},{"RecvDepth",encode(v.Rz,v.is_Linspace_Rz,v.NRz,"NRz")},{"RecvAzim",encode(v.Ro,v.is_Linspace_Ro,v.NRo,"NRo")}}; }
void from_json(const OpenOcean_json &in, Position &v) { const auto decode=[](const OpenOcean_json &item,const char *countName,const char *field,Eigen::VectorXd &values,int &count,bool &linspace){ if(item.is_object()&&item.contains("start")){ count=item.at(countName).get<int>(); if(count<1) throw std::runtime_error(std::string(field)+" count must be positive"); values=Eigen::VectorXd::LinSpaced(count,item.at("start").get<double>(),item.at("end").get<double>()); linspace=true; } else { values=readRealVector(item,field); count=static_cast<int>(values.size()); linspace=false; } }; if(in.contains("SrcDepth")){ decode(in.at("SrcDepth"),"NSz","SrcDepth",v.Sz,v.NSz,v.is_Linspace_Sz); decode(in.at("RecvRange"),"NRr","RecvRange",v.Rr,v.NRr,v.is_Linspace_Rr); decode(in.at("RecvDepth"),"NRz","RecvDepth",v.Rz,v.NRz,v.is_Linspace_Rz); decode(in.at("RecvAzim"),"NRo","RecvAzim",v.Ro,v.NRo,v.is_Linspace_Ro); } else { v.Sz=readRealVector(in.at("Sz"),"Sz"); v.Rr=readRealVector(in.at("Rr"),"Rr"); v.Rz=readRealVector(in.at("Rz"),"Rz"); v.Ro=readRealVector(in.at("Ro"),"Ro"); v.NSz=in.value("NSz",static_cast<int>(v.Sz.size())); v.NRr=in.value("NRr",static_cast<int>(v.Rr.size())); v.NRz=in.value("NRz",static_cast<int>(v.Rz.size())); v.NRo=in.value("NRo",static_cast<int>(v.Ro.size())); v.is_Linspace_Rr=in.value("is_Linspace_Rr",false); v.is_Linspace_Rz=in.value("is_Linspace_Rz",false); v.is_Linspace_Sz=in.value("is_Linspace_Sz",false); v.is_Linspace_Ro=in.value("is_Linspace_Ro",false); } in.at("GridType").get_to(v.GridType); v.NRz_per_range=in.value("NRz_per_range",v.GridType==Grid_Mode::MODE_I_Irregular?1:v.NRz); v.Delta_r=in.value("Delta_r",0.0); checkPosition(v); }
void to_json(OpenOcean_json &out, const ReflectionCoef &v) { out={{"theta",v.theta},{"R",v.R},{"phi",v.phi}}; }
void from_json(const OpenOcean_json &in, ReflectionCoef &v) { v.theta=in.at("theta").get<double>(); v.R=in.at("R").get<double>(); v.phi=in.at("phi").get<double>(); if(!std::isfinite(v.theta)||!std::isfinite(v.R)||!std::isfinite(v.phi)||v.R<0.0) throw std::runtime_error("reflection coefficient values are invalid"); }
void to_json(OpenOcean_json &out, const InternalReflectionCoefInfo &v) { if(!v.isSet){ out={{"isSet",false}}; return; } OpenOcean_json fReal=OpenOcean_json::array(),fImag=OpenOcean_json::array(),gReal=OpenOcean_json::array(),gImag=OpenOcean_json::array(); for(Eigen::Index i=0;i<v.fTab.size();++i){ fReal.push_back(v.fTab[i].real()); fImag.push_back(v.fTab[i].imag()); gReal.push_back(v.gTab[i].real()); gImag.push_back(v.gTab[i].imag()); } out={{"isSet",true},{"freq",v.freq},{"xTab",realVector(v.xTab)},{"fReal",fReal},{"fImag",fImag},{"gReal",gReal},{"gImag",gImag},{"iTab",intVector(v.iTab)}}; }
void from_json(const OpenOcean_json &in, InternalReflectionCoefInfo &v) { if(!in.is_object()||!in.value("isSet",false)){ v=InternalReflectionCoefInfo{}; return; } v.freq=in.at("freq").get<double>(); v.xTab=readRealVector(in.at("xTab"),"xTab"); v.iTab=readIntVector(in.at("iTab"),"iTab"); if(in.contains("fTab")){ v.fTab=readComplexVector(in.at("fTab"),"fTab"); v.gTab=readComplexVector(in.at("gTab"),"gTab"); } else { const Eigen::VectorXd fr=readRealVector(in.at("fReal"),"fReal"),fi=readRealVector(in.at("fImag"),"fImag"),gr=readRealVector(in.at("gReal"),"gReal"),gi=readRealVector(in.at("gImag"),"gImag"); if(fr.size()!=fi.size()||fr.size()!=gr.size()||fr.size()!=gi.size()) throw std::runtime_error("internal reflection component dimensions mismatch"); v.fTab.resize(fr.size()); v.gTab.resize(fr.size()); for(Eigen::Index i=0;i<fr.size();++i){ v.fTab[i]={fr[i],fi[i]}; v.gTab[i]={gr[i],gi[i]}; } } v.isSet=true; if(v.xTab.size()!=v.fTab.size()||v.xTab.size()!=v.gTab.size()||v.xTab.size()!=v.iTab.size()) throw std::runtime_error("internal reflection dimensions mismatch"); }
void to_json(OpenOcean_json &out, const ReflectionCoefInfo &v) { const auto encode=[&](const Eigen::Matrix<ReflectionCoef,1,Eigen::Dynamic> &values){ OpenOcean_json r=OpenOcean_json::array(),phi=OpenOcean_json::array(),theta=OpenOcean_json::array(); const double factor=v.isDeg?180.0/pi:1.0; for(Eigen::Index i=0;i<values.size();++i){ r.push_back(values[i].R); phi.push_back(values[i].phi*factor); theta.push_back(values[i].theta); } return OpenOcean_json{{"R",r},{"phi",phi},{"theta",theta}}; }; out={{"RTop",encode(v.RTop)},{"RBot",encode(v.RBot)},{"IRC",v.IRC}}; }
void from_json(const OpenOcean_json &in, ReflectionCoefInfo &v) { const auto decode=[](const OpenOcean_json &item){ Eigen::Matrix<ReflectionCoef,1,Eigen::Dynamic> result; if(item.is_array()){ result.resize(1,static_cast<Eigen::Index>(item.size())); for(Eigen::Index i=0;i<result.size();++i) item.at(static_cast<std::size_t>(i)).get_to(result[i]); return result; } const Eigen::VectorXd r=readRealVector(item.at("R"),"reflection R"),phi=readRealVector(item.at("phi"),"reflection phi"),theta=readRealVector(item.at("theta"),"reflection theta"); if(r.size()!=phi.size()||r.size()!=theta.size()) throw std::runtime_error("reflection coefficient dimensions mismatch"); result.resize(1,r.size()); for(Eigen::Index i=0;i<r.size();++i) { if(!std::isfinite(theta[i])||!std::isfinite(r[i])||!std::isfinite(phi[i])||r[i]<0.0||(i>0&&theta[i]<=theta[i-1])) throw std::runtime_error("reflection coefficients must be finite and angle ordered"); result[i]={theta[i],r[i],phi[i]}; } return result; }; v.RTop=decode(in.at("RTop")); v.RBot=decode(in.at("RBot")); if(in.contains("IRC")) in.at("IRC").get_to(v.IRC); else v.IRC=InternalReflectionCoefInfo{}; v.isDeg=in.value("isDeg",false); }
void to_json(OpenOcean_json &out, const SrcBmPat &v) { if(!v.isSet){ out=nullptr; return; } out={{"NSBPPts",v.NSBPPts},{"theta",realVector(v.theta)},{"pat",realVector(v.pat)}}; }
void from_json(const OpenOcean_json &in, SrcBmPat &v) { if(in.is_null()||!in.is_object()||!in.contains("NSBPPts")||(in.contains("isSet")&&!in.at("isSet").get<bool>())){ v=SrcBmPat{}; return; } v.theta=readRealVector(in.at("theta"),"SBP theta"); v.pat=readRealVector(in.at("pat"),"SBP pat"); v.NSBPPts=in.value("NSBPPts",static_cast<int>(v.theta.size())); v.isSet=true; if(v.NSBPPts!=v.theta.size()||v.theta.size()!=v.pat.size()||v.theta.size()<2||!v.theta.allFinite()||!v.pat.allFinite()||(v.theta.tail(v.theta.size()-1)-v.theta.head(v.theta.size()-1)).minCoeff()<=0.0) throw std::runtime_error("SBP must be finite, dimensionally consistent, and strictly ordered"); }
void to_json(OpenOcean_json &out, const BdryType &v) { out={{"Top",v.Top},{"Bot",v.Bot}}; }
void from_json(const OpenOcean_json &in, BdryType &v) { in.at("Top").get_to(v.Top); in.at("Bot").get_to(v.Bot); }
void to_json(OpenOcean_json &out, const FreqInfo &v) { out["Nfreq"]=v.Nfreq; out["freq"]=v.Nfreq>1?realVector(v.freqvec):OpenOcean_json(v.freq); }
void from_json(const OpenOcean_json &in, FreqInfo &v) { v.Nfreq=in.value("Nfreq",1); if(in.at("freq").is_array()){ v.freqvec=readRealVector(in.at("freq"),"freq"); if(v.Nfreq!=v.freqvec.size()) throw std::runtime_error("Nfreq does not match frequency vector"); v.freq=v.freqvec.size()?v.freqvec[0]:0.0; } else { v.freq=in.at("freq").get<double>(); if(in.contains("freqvec")) v.freqvec=readRealVector(in.at("freqvec"),"freqvec"); else v.freqvec=Eigen::VectorXd::Constant(1,v.freq); if(v.Nfreq!=1) throw std::runtime_error("multi-frequency input requires a frequency array"); } if(v.Nfreq<1||v.freqvec.size()!=v.Nfreq||!v.freqvec.allFinite()||(v.freqvec.array()<=0.0).any()||!std::isfinite(v.freq)||v.freq<=0.0) throw std::runtime_error("frequencies must be finite, positive, and dimensionally consistent"); }
void to_json(OpenOcean_json &out, const MeshParams &v) { out={{"NSets",v.NSets},{"NV",intVector(v.NV)}}; }
void from_json(const OpenOcean_json &in, MeshParams &v) { v.NSets=in.at("NSets").get<int>(); v.NV=readIntVector(in.at("NV"),"NV"); if(v.NSets<1||v.NV.size()!=v.NSets) throw std::runtime_error("mesh dimensions mismatch"); }

namespace ssp
{
void to_json(OpenOcean_json &out, const SSPLayer &v) { out={{"npts",v.npts},{"nmesh",v.nmesh},{"beta",v.beta},{"ft",v.ft},{"sigma",v.sigma},{"z",realVector(v.z)},{"rho",realVector(v.rho)},{"alphaR",realVector(v.alphaR)},{"alphaI",realVector(v.alphaI)},{"betaR",realVector(v.betaR)},{"betaI",realVector(v.betaI)},{"Material",v.Material}}; }
void from_json(const OpenOcean_json &in, SSPLayer &v) { v.npts=in.at("npts").get<int>(); v.nmesh=in.at("nmesh").get<int>(); v.beta=in.at("beta").get<double>(); v.ft=in.at("ft").get<double>(); v.sigma=in.at("sigma").get<double>(); v.z=readRealVector(in.at("z"),"z"); v.rho=readRealVector(in.at("rho"),"rho"); v.alphaR=readRealVector(in.at("alphaR"),"alphaR"); v.alphaI=readRealVector(in.at("alphaI"),"alphaI"); v.betaR=readRealVector(in.at("betaR"),"betaR"); v.betaI=readRealVector(in.at("betaI"),"betaI"); in.at("Material").get_to(v.Material); if(!v.is_valid()||v.nmesh<1||!std::isfinite(v.beta)||!std::isfinite(v.ft)||!std::isfinite(v.sigma)||!v.z.allFinite()||!v.rho.allFinite()||!v.alphaR.allFinite()||!v.alphaI.allFinite()||!v.betaR.allFinite()||!v.betaI.allFinite()||(v.z.tail(v.z.size()-1)-v.z.head(v.z.size()-1)).minCoeff()<=0.0||(v.rho.array()<=0.0).any()||(v.alphaR.array()<=0.0).any()||(v.betaR.array()<0.0).any()||(v.alphaI.array()<0.0).any()||(v.betaI.array()<0.0).any()) throw std::runtime_error("invalid SSP layer dimensions or values"); }
void to_json(OpenOcean_json &out, const Range_Independent_Area &v) { out={{"Title",v.Title},{"SSPType",v.SSPType},{"Range",v.Range},{"enableRootRestarts",v.enableRootRestarts},{"layers",v.layers},{"HSTop",v.HSTop},{"HSBot",v.HSBot}}; }
void from_json(const OpenOcean_json &in, Range_Independent_Area &v) { v.Title=in.value("Title",""); in.at("SSPType").get_to(v.SSPType); v.Range=in.at("Range").get<double>(); v.enableRootRestarts=in.value("enableRootRestarts",false); v.layers.clear(); for(const auto &item:in.at("layers")) v.layers.push_back(item.get<SSPLayer>()); in.at("HSTop").get_to(v.HSTop); in.at("HSBot").get_to(v.HSBot); }
void to_json(OpenOcean_json &out, const FlattenedData &v) { out={{"NPts",intVector(v.NPts)},{"NMesh",intVector(v.NMesh)},{"beta",realVector(v.beta)},{"ft",realVector(v.ft)},{"sigma",realVector(v.sigma)},{"z",realVector(v.z)},{"rho",realVector(v.rho)},{"alphaR",realVector(v.alphaR)},{"alphaI",realVector(v.alphaI)},{"betaR",realVector(v.betaR)},{"betaI",realVector(v.betaI)}}; }
void from_json(const OpenOcean_json &in, FlattenedData &v) { v.NPts=readIntVector(in.at("NPts"),"NPts"); v.NMesh=readIntVector(in.at("NMesh"),"NMesh"); v.beta=readRealVector(in.at("beta"),"beta"); v.ft=readRealVector(in.at("ft"),"ft"); v.sigma=readRealVector(in.at("sigma"),"sigma"); v.z=readRealVector(in.at("z"),"z"); v.rho=readRealVector(in.at("rho"),"rho"); v.alphaR=readRealVector(in.at("alphaR"),"alphaR"); v.alphaI=readRealVector(in.at("alphaI"),"alphaI"); v.betaR=readRealVector(in.at("betaR"),"betaR"); v.betaI=readRealVector(in.at("betaI"),"betaI"); }
void to_json(OpenOcean_json &out, const SSPStructure &v) { out={{"NMedia",v.NMedia},{"FirstAcoustic",v.FirstAcoustic},{"LastAcoustic",v.LastAcoustic},{"beta",realVector(v.beta)},{"ft",realVector(v.ft)},{"sigma",realVector(v.sigma)},{"Material",v.Material},{"SSPType",v.SSPType},{"depth",realVector(v.depth)},{"NPts",intVector(v.NPts)},{"offset",intVector(v.offset)},{"interp_offset",intVector(v.interp_offset)},{"NMesh",intVector(v.NMesh)},{"z",realVector(v.z)},{"alphaR",realVector(v.alphaR)},{"alphaI",realVector(v.alphaI)},{"betaR",realVector(v.betaR)},{"betaI",realVector(v.betaI)},{"rho",realVector(v.rho)},{"cp",complexVector(v.cp)},{"cs",complexVector(v.cs)},{"row_size",v.row_size},{"cspline",complexMatrix(v.cspline)},{"cCoef",complexMatrix(v.cCoef)},{"csWork",complexMatrix(v.csWork)},{"cpCoef",complexMatrix(v.cpCoef)},{"csCoef",complexMatrix(v.csCoef)},{"rhoCoef",complexMatrix(v.rhoCoef)},{"cpSpline",complexMatrix(v.cpSpline)},{"csSpline",complexMatrix(v.csSpline)},{"rhoSpline",complexMatrix(v.rhoSpline)},{"HSTop",v.HSTop},{"HSBot",v.HSBot}}; }
void from_json(const OpenOcean_json &in, SSPStructure &v) { v.NMedia=in.at("NMedia").get<int>(); v.FirstAcoustic=in.at("FirstAcoustic").get<int>(); v.LastAcoustic=in.at("LastAcoustic").get<int>(); v.beta=readRealVector(in.at("beta"),"beta"); v.ft=readRealVector(in.at("ft"),"ft"); v.sigma=readRealVector(in.at("sigma"),"sigma"); v.Material.clear(); for(const auto &item:in.at("Material")) v.Material.push_back(item.get<Media_Mode>()); in.at("SSPType").get_to(v.SSPType); v.depth=readRealVector(in.at("depth"),"depth"); v.NPts=readIntVector(in.at("NPts"),"NPts"); v.offset=readIntVector(in.at("offset"),"offset"); v.interp_offset=readIntVector(in.at("interp_offset"),"interp_offset"); v.NMesh=readIntVector(in.at("NMesh"),"NMesh"); v.z=readRealVector(in.at("z"),"z"); v.alphaR=readRealVector(in.at("alphaR"),"alphaR"); v.alphaI=readRealVector(in.at("alphaI"),"alphaI"); v.betaR=readRealVector(in.at("betaR"),"betaR"); v.betaI=readRealVector(in.at("betaI"),"betaI"); v.rho=readRealVector(in.at("rho"),"rho"); v.cp=readComplexVector(in.at("cp"),"cp"); v.cs=readComplexVector(in.at("cs"),"cs"); v.row_size=in.at("row_size").get<int>(); v.cspline=readComplexMatrix(in.at("cspline"),"cspline"); v.cCoef=readComplexMatrix(in.at("cCoef"),"cCoef"); v.csWork=readComplexMatrix(in.at("csWork"),"csWork"); v.cpCoef=readComplexMatrix(in.at("cpCoef"),"cpCoef"); v.csCoef=readComplexMatrix(in.at("csCoef"),"csCoef"); v.rhoCoef=readComplexMatrix(in.at("rhoCoef"),"rhoCoef"); v.cpSpline=readComplexMatrix(in.at("cpSpline"),"cpSpline"); v.csSpline=readComplexMatrix(in.at("csSpline"),"csSpline"); v.rhoSpline=readComplexMatrix(in.at("rhoSpline"),"rhoSpline"); in.at("HSTop").get_to(v.HSTop); in.at("HSBot").get_to(v.HSBot); }
}

void to_json(OpenOcean_json &out, const OOKC_parameters &p)
{
    const int profileCount = p.sspInput.empty() ? 0 : p.NProf;
    const Eigen::VectorXd profileRanges = p.sspInput.empty() ? Eigen::VectorXd{} : p.RProf;
    out={{"schema","OpenOcean-Krakenc.parameters"},{"schema_version",1},{"schemaVersion",1},{"Title",p.Title},{"freqinfo",p.freqinfo},{"AttenUnit",p.AttenUnit},{"Pos",p.Pos},{"MLimit",p.MLimit},{"NProf",profileCount},{"RProf",realVector(profileRanges)},{"hasModePos",p.hasModePos},{"sspInput",p.sspInput},{"ReflectionCoef",p.ReflectionCoef},{"SBP",p.SBP},{"is_Velocity",p.is_Velocity},{"cLow",p.cLow},{"cHigh",p.cHigh},{"Rmax",p.Rmax},{"SourceType",p.SourceType},{"mesh",p.mesh},{"RunMode",p.runMode},{"CoherenceType",p.coherenceType},{"ModeType",p.modeType}};
    if(p.hasModePos) out["ModePos"]=p.ModePos;
    out["paths"]={{"env",p.envPath},{"flp",p.flpPath},{"mod",p.modPath},{"shd",p.shdPath}};
}

void from_json(const OpenOcean_json &in, OOKC_parameters &p)
{
    if(in.contains("schema")&&in.at("schema").get<std::string>()!="OpenOcean-Krakenc.parameters") throw std::runtime_error("unsupported JSON schema");
    const int version=in.contains("schema_version")?in.at("schema_version").get<int>():in.value("schemaVersion",1);
    if(in.contains("schema_version")&&in.contains("schemaVersion")&&in.at("schemaVersion").get<int>()!=version) throw std::runtime_error("conflicting JSON schema versions");
    if(version!=1) throw std::runtime_error("unsupported JSON schema version");
    OOKC_parameters parsed;
    parsed.Title=in.at("Title").get<std::string>(); in.at("freqinfo").get_to(parsed.freqinfo); in.at("AttenUnit").get_to(parsed.AttenUnit); in.at("Pos").get_to(parsed.Pos);
    parsed.MLimit=in.value("MLimit",9999); parsed.hasModePos=in.value("hasModePos",false); parsed.ModePos=(parsed.hasModePos&&in.contains("ModePos"))?in.at("ModePos").get<Position>():parsed.Pos;
    for(const auto &item:in.at("sspInput")) parsed.sspInput.push_back(item.get<ssp::Range_Independent_Area>());
    parsed.NProf=in.value("NProf",static_cast<int>(parsed.sspInput.size()));
    if(in.contains("RProf")) parsed.RProf=readRealVector(in.at("RProf"),"RProf"); else { parsed.RProf.resize(parsed.NProf); for(Eigen::Index i=0;i<parsed.RProf.size();++i) parsed.RProf[i]=parsed.sspInput.at(static_cast<std::size_t>(i)).Range; }
    if(parsed.NProf!=static_cast<int>(parsed.sspInput.size())||parsed.RProf.size()!=parsed.NProf) throw std::runtime_error("NProf/RProf must match SSP profiles");
    for(Eigen::Index i=0;i<parsed.RProf.size();++i) { if(!std::isfinite(parsed.RProf[i])||(i==0&&std::abs(parsed.RProf[i])>1e-9)||(i>0&&parsed.RProf[i]<=parsed.RProf[i-1])) throw std::runtime_error("RProf must start at zero and increase"); parsed.sspInput[static_cast<std::size_t>(i)].Range=parsed.RProf[i]; }
    in.at("ReflectionCoef").get_to(parsed.ReflectionCoef); in.at("SBP").get_to(parsed.SBP); parsed.is_Velocity=in.at("is_Velocity").get<bool>(); parsed.cLow=in.at("cLow").get<double>(); parsed.cHigh=in.at("cHigh").get<double>(); parsed.Rmax=in.at("Rmax").get<double>(); in.at("SourceType").get_to(parsed.SourceType); if(in.contains("mesh")) in.at("mesh").get_to(parsed.mesh); in.at("RunMode").get_to(parsed.runMode); in.at("CoherenceType").get_to(parsed.coherenceType); in.at("ModeType").get_to(parsed.modeType);
    if(parsed.MLimit<1) throw std::runtime_error("MLimit must be positive");
    if(!std::isfinite(parsed.Rmax)||parsed.Rmax<0.0) throw std::runtime_error("Rmax must be finite and nonnegative");
    if((!parsed.sspInput.empty()||parsed.cLow!=0.0||parsed.cHigh!=0.0)&&(!std::isfinite(parsed.cLow)||!std::isfinite(parsed.cHigh)||parsed.cLow<=0.0||parsed.cHigh<=parsed.cLow)) throw std::runtime_error("phase speed interval is invalid");
    for(const auto &area:parsed.sspInput) for(const auto &layer:area.layers) { const bool shear=(layer.betaR.array().abs()>1.0e-12).any(); if(shear!=(layer.Material==Media_Mode::MODE_E_Elastic)) throw std::runtime_error("SSP Material does not match shear speed"); }
    if(in.contains("paths")) { const auto &paths=in.at("paths"); parsed.envPath=paths.value("env",""); parsed.flpPath=paths.value("flp",""); parsed.modPath=paths.value("mod",""); parsed.shdPath=paths.value("shd",""); }
    if(!parsed.sspInput.empty())
    {
        OOKC_parameters semanticCandidate=parsed;
        semanticCandidate.freqinfo.Nfreq=1;
        semanticCandidate.freqinfo.freq=parsed.freqinfo.freqvec[0];
        semanticCandidate.freqinfo.freqvec=Eigen::VectorXd::Constant(
            1,semanticCandidate.freqinfo.freq);
        validatePublicParameterSemantics(semanticCandidate,
                                         semanticCandidate.runMode);
    }
    p=std::move(parsed);
}

bool read_json_file(const std::string &path, OOKC_parameters &params)
{
    try { std::ifstream stream(path); if(!stream) return false; OpenOcean_json document; stream>>document; OOKC_parameters candidate; from_json(document,candidate); params=std::move(candidate); return true; } catch(...) { return false; }
}
bool write_json_file(const std::string &path, const OOKC_parameters &params)
{
    try { OpenOcean_json document; to_json(document,params); std::ofstream stream(path); if(!stream) return false; stream<<document.dump(4)<<'\n'; return static_cast<bool>(stream); } catch(...) { return false; }
}
std::string parameters_to_json_string(const OOKC_parameters &params) { OpenOcean_json document; to_json(document,params); return document.dump(4); }
}
