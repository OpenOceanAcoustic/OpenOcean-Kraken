#include "module/ParameterAdapters.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

namespace OpenOceanKrakenc
{
namespace
{
char attenuationCode(AttenuationUnit unit)
{
    switch (unit)
    {
    case AttenuationUnit::MODE_F_dB_per_m_kHz: return 'F';
    case AttenuationUnit::MODE_L_params_lose: return 'L';
    case AttenuationUnit::MODE_M_dB_per_m: return 'M';
    case AttenuationUnit::MODE_m_dB_per_m: return 'm';
    case AttenuationUnit::MODE_N_Nepers_per_m: return 'N';
    case AttenuationUnit::MODE_Q_Quality_Factor: return 'Q';
    case AttenuationUnit::MODE_W_db_per_lambda: return 'W';
    }
    throw std::invalid_argument("unsupported attenuation unit");
}

AttenuationUnit attenuationUnit(char code)
{
    switch (code)
    {
    case 'F': return AttenuationUnit::MODE_F_dB_per_m_kHz;
    case 'L': return AttenuationUnit::MODE_L_params_lose;
    case 'M': return AttenuationUnit::MODE_M_dB_per_m;
    case 'm': return AttenuationUnit::MODE_m_dB_per_m;
    case 'N': return AttenuationUnit::MODE_N_Nepers_per_m;
    case 'Q': return AttenuationUnit::MODE_Q_Quality_Factor;
    case 'W': return AttenuationUnit::MODE_W_db_per_lambda;
    default: throw std::invalid_argument("unsupported ENV attenuation unit");
    }
}

AcousticInterpolation interpolation(SSP_Mode type)
{
    switch (type)
    {
    case SSP_Mode::MODE_N_n2Linear: return AcousticInterpolation::N2Linear;
    case SSP_Mode::MODE_C_cLinear: return AcousticInterpolation::CLinear;
    case SSP_Mode::MODE_P_cPCHIP: return AcousticInterpolation::Pchip;
    case SSP_Mode::MODE_S_cCubic: return AcousticInterpolation::CubicSpline;
    case SSP_Mode::MODE_A_Analytic:
        throw std::invalid_argument(
            "analytic SSP interpolation is not implemented");
    }
    throw std::invalid_argument("unsupported public SSP interpolation");
}

SSP_Mode interpolation(AcousticInterpolation type)
{
    switch (type)
    {
    case AcousticInterpolation::N2Linear:
        return SSP_Mode::MODE_N_n2Linear;
    case AcousticInterpolation::CLinear:
        return SSP_Mode::MODE_C_cLinear;
    case AcousticInterpolation::Pchip:
        return SSP_Mode::MODE_P_cPCHIP;
    case AcousticInterpolation::CubicSpline:
        return SSP_Mode::MODE_S_cCubic;
    }
    throw std::invalid_argument("unsupported internal SSP interpolation");
}

AcousticBoundaryType boundaryType(BC_Mode type)
{
    switch (type)
    {
    case BC_Mode::MODE_V_Vacuum: return AcousticBoundaryType::Vacuum;
    case BC_Mode::MODE_R_Rigid: return AcousticBoundaryType::Rigid;
    case BC_Mode::MODE_A_Half_space: return AcousticBoundaryType::HalfSpace;
    case BC_Mode::MODE_F_File: return AcousticBoundaryType::ReflectionCoefficient;
    case BC_Mode::MODE_P_Precomputed: return AcousticBoundaryType::InternalReflection;
    case BC_Mode::MODE_G_Grain:
        throw std::invalid_argument("Krakenc grain boundary is not implemented");
    }
    throw std::invalid_argument("unsupported boundary type");
}

BC_Mode boundaryType(AcousticBoundaryType type)
{
    switch (type)
    {
    case AcousticBoundaryType::Vacuum: return BC_Mode::MODE_V_Vacuum;
    case AcousticBoundaryType::Rigid: return BC_Mode::MODE_R_Rigid;
    case AcousticBoundaryType::HalfSpace: return BC_Mode::MODE_A_Half_space;
    case AcousticBoundaryType::ReflectionCoefficient: return BC_Mode::MODE_F_File;
    case AcousticBoundaryType::InternalReflection: return BC_Mode::MODE_P_Precomputed;
    }
    throw std::invalid_argument("unsupported internal boundary type");
}

AcousticBoundary makeBoundary(
    const HSInfo &source,
    const Eigen::Matrix<ReflectionCoef, 1, Eigen::Dynamic> &reflection,
    const InternalReflectionCoefInfo *internal,
    bool reflectionIsRadians)
{
    AcousticBoundary result;
    result.type = boundaryType(source.BC);
    result.depth = source.Depth;
    result.cp = source.alphaR;
    result.cs = source.betaR;
    result.rho = source.rho;
    result.alphaP = source.alphaI;
    result.alphaS = source.betaI;
    result.roughnessRms = source.sigma;
    result.attenuationPower = source.beta;
    result.transitionFrequency = source.ft;
    for (Eigen::Index i = 0; i < reflection.size(); ++i)
    {
        const double phaseRadians = reflectionIsRadians
                                        ? reflection[i].phi
                                        : reflection[i].phi * std::acos(-1.0) / 180.0;
        result.reflectionSamples.push_back(
            {reflection[i].theta, reflection[i].R, phaseRadians});
    }
    if (internal && internal->isSet)
    {
        const Eigen::Index count = internal->xTab.size();
        if (internal->fTab.size() != count || internal->gTab.size() != count ||
            internal->iTab.size() != count)
        {
            throw std::invalid_argument("internal reflection table dimensions mismatch");
        }
        for (Eigen::Index i = 0; i < count; ++i)
        {
            result.internalSamples.push_back(
                {internal->xTab[i], internal->fTab[i], internal->gTab[i], internal->iTab[i]});
        }
    }
    return result;
}

HSInfo makeHalfSpace(const AcousticBoundary &source)
{
    HSInfo result;
    result.BC = boundaryType(source.type);
    result.Depth = source.depth;
    result.alphaR = source.cp;
    result.alphaI = source.alphaP;
    result.betaR = source.cs;
    result.betaI = source.alphaS;
    result.sigma = source.roughnessRms;
    result.beta = source.attenuationPower;
    result.ft = source.transitionFrequency;
    result.cp = {source.cp, source.alphaP};
    result.cs = {source.cs, source.alphaS};
    result.rho = source.rho;
    return result;
}

template <typename Scalar>
std::vector<double> toVector(const Eigen::Matrix<Scalar, Eigen::Dynamic, 1> &values)
{
    return std::vector<double>(values.data(), values.data() + values.size());
}

Eigen::VectorXd toEigen(const std::vector<double> &values)
{
    Eigen::VectorXd result(static_cast<Eigen::Index>(values.size()));
    for (Eigen::Index i = 0; i < result.size(); ++i)
    {
        result[i] = values[static_cast<std::size_t>(i)];
    }
    return result;
}

char sourceCode(Source_Mode mode)
{
    switch (mode)
    {
    case Source_Mode::MODE_R_Point: return 'R';
    case Source_Mode::MODE_X_Line: return 'X';
    case Source_Mode::MODE_S_ScaledCylindrical: return 'S';
    }
    throw std::invalid_argument("unsupported source type");
}

Source_Mode sourceMode(char code)
{
    switch (code)
    {
    case 'R': return Source_Mode::MODE_R_Point;
    case 'X': return Source_Mode::MODE_X_Line;
    case 'S': return Source_Mode::MODE_S_ScaledCylindrical;
    default: throw std::invalid_argument("unsupported FLP source type");
    }
}

void requireNonempty(const Eigen::VectorXd &values, const char *name)
{
    if (values.size() < 1)
    {
        throw std::invalid_argument(std::string(name) + " must not be empty");
    }
}
}

void validatePublicParameterSemantics(const OOKC_parameters &params,
                                      Run_Mode requestedRun)
{
    if (params.freqinfo.freqvec.size() > 0)
    {
        if ((params.freqinfo.Nfreq > 0 &&
             params.freqinfo.Nfreq != params.freqinfo.freqvec.size()) ||
            !params.freqinfo.freqvec.allFinite() ||
            (params.freqinfo.freqvec.array() <= 0.0).any())
        {
            throw std::invalid_argument(
                "frequencies must be finite, positive, and dimensionally consistent");
        }
    }
    const double frequency = params.freqinfo.freqvec.size() == 1
                                 ? params.freqinfo.freqvec[0]
                                 : params.freqinfo.freq;
    if (!std::isfinite(frequency) || !(frequency > 0.0))
    {
        throw std::invalid_argument("frequency must be positive");
    }
    if (!std::isfinite(params.Rmax) || params.Rmax < 0.0)
    {
        throw std::invalid_argument("Rmax must be finite and nonnegative");
    }
    if (params.sspInput.empty())
    {
        throw std::invalid_argument("SSP input must not be empty");
    }
    if (params.NProf != 0 && params.NProf != static_cast<int>(params.sspInput.size()))
    {
        throw std::invalid_argument("NProf does not match SSP profile count");
    }
    if (params.RProf.size() != 0 && params.RProf.size() != static_cast<Eigen::Index>(params.sspInput.size()))
    {
        throw std::invalid_argument("RProf does not match SSP profile count");
    }
    if (params.RProf.size() != 0)
    {
        for (Eigen::Index profile = 0; profile < params.RProf.size(); ++profile)
        {
            const double rangeKm = params.RProf[profile];
            if (!std::isfinite(rangeKm) ||
                (profile == 0 && std::abs(rangeKm) > 1.0e-9) ||
                (profile > 0 && rangeKm <= params.RProf[profile - 1]))
            {
                throw std::invalid_argument(
                    "RProf must be finite, start at zero, and strictly increase in km");
            }
        }
    }
    for (const ssp::Range_Independent_Area &area : params.sspInput)
    {
        const auto validateBoundary = [](const HSInfo &boundary,
                                         const char *name) {
            const bool finite = std::isfinite(boundary.alphaR) &&
                                std::isfinite(boundary.alphaI) &&
                                std::isfinite(boundary.betaR) &&
                                std::isfinite(boundary.betaI) &&
                                std::isfinite(boundary.beta) &&
                                std::isfinite(boundary.ft) &&
                                std::isfinite(boundary.rho) &&
                                std::isfinite(boundary.Depth) &&
                                std::isfinite(boundary.cp.real()) &&
                                std::isfinite(boundary.cp.imag()) &&
                                std::isfinite(boundary.cs.real()) &&
                                std::isfinite(boundary.cs.imag());
            if (!finite)
            {
                throw std::invalid_argument(std::string(name) +
                                            " boundary values must be finite");
            }
            if (boundary.BC == BC_Mode::MODE_A_Half_space &&
                (!(boundary.alphaR > 0.0) || boundary.betaR < 0.0 ||
                 !(boundary.rho > 0.0) || boundary.alphaI < 0.0 ||
                 boundary.betaI < 0.0))
            {
                throw std::invalid_argument(std::string(name) +
                                            " halfspace properties are invalid");
            }
        };
        validateBoundary(area.HSTop, "top");
        validateBoundary(area.HSBot, "bottom");
        const auto validateReflection = [](const auto &samples,
                                           const char *name) {
            if (samples.size() < 2)
            {
                throw std::invalid_argument(std::string(name) +
                                            " reflection table requires two samples");
            }
            for (Eigen::Index i = 0; i < samples.size(); ++i)
            {
                if (!std::isfinite(samples[i].theta) ||
                    !std::isfinite(samples[i].R) ||
                    !std::isfinite(samples[i].phi) || samples[i].R < 0.0 ||
                    (i > 0 && samples[i].theta <= samples[i - 1].theta))
                {
                    throw std::invalid_argument(std::string(name) +
                                                " reflection table is invalid");
                }
            }
        };
        if (area.HSTop.BC == BC_Mode::MODE_F_File)
        {
            validateReflection(params.ReflectionCoef.RTop, "top");
        }
        if (area.HSBot.BC == BC_Mode::MODE_F_File)
        {
            validateReflection(params.ReflectionCoef.RBot, "bottom");
        }
        if (area.layers.empty())
        {
            throw std::invalid_argument("each SSP profile must contain a layer");
        }
        for (const ssp::SSPLayer &layer : area.layers)
        {
            if (!layer.is_valid() || layer.nmesh < 1)
            {
                throw std::invalid_argument("SSP layer dimensions or mesh are invalid");
            }
            if (!std::isfinite(layer.beta) || !std::isfinite(layer.ft) ||
                !std::isfinite(layer.sigma) || !layer.z.allFinite() ||
                !layer.rho.allFinite() || !layer.alphaR.allFinite() ||
                !layer.alphaI.allFinite() || !layer.betaR.allFinite() ||
                !layer.betaI.allFinite())
            {
                throw std::invalid_argument("SSP layer values must be finite");
            }
            if ((layer.z.tail(layer.z.size() - 1) -
                 layer.z.head(layer.z.size() - 1)).minCoeff() <= 0.0)
            {
                throw std::invalid_argument("SSP depths must strictly increase");
            }
            if ((layer.alphaR.array() <= 0.0).any() ||
                (layer.rho.array() <= 0.0).any() ||
                (layer.betaR.array() < 0.0).any() ||
                (layer.alphaI.array() < 0.0).any() ||
                (layer.betaI.array() < 0.0).any())
            {
                throw std::invalid_argument(
                    "SSP sound speeds, density, or attenuation are invalid");
            }
            const bool allCsZero =
                (layer.betaR.array() == 0.0).all();
            const bool allCsPositive =
                (layer.betaR.array() > 0.0).all();
            if (!allCsZero && !allCsPositive)
            {
                throw std::invalid_argument(
                    "SSP shear sound speeds must be either all zero or all positive");
            }
            const bool materialMatches =
                (allCsZero &&
                 layer.Material == Media_Mode::MODE_A_Acoustic) ||
                (allCsPositive &&
                 layer.Material == Media_Mode::MODE_E_Elastic);
            if (!materialMatches)
            {
                throw std::invalid_argument(
                    "SSP layer Material must agree with its shear sound speed");
            }
        }
    }
    if (requestedRun != Run_Mode::MODE_F_Field &&
        (!std::isfinite(params.cLow) || !std::isfinite(params.cHigh) ||
         !(params.cLow > 0.0) || !(params.cHigh > params.cLow)))
    {
        throw std::invalid_argument("phase speed interval is invalid");
    }
    const Position &modePosition = params.hasModePos ? params.ModePos : params.Pos;
    requireNonempty(modePosition.Sz, "source depths");
    requireNonempty(modePosition.Rz, "receiver depths");
    if (!modePosition.Sz.allFinite() || !modePosition.Rz.allFinite())
    {
        throw std::invalid_argument("mode source and receiver depths must be finite");
    }
    if (requestedRun != Run_Mode::MODE_M_Modes)
    {
        requireNonempty(params.Pos.Sz, "field source depths");
        requireNonempty(params.Pos.Rz, "field receiver depths");
        requireNonempty(params.Pos.Rr, "receiver ranges");
        if (!params.Pos.Sz.allFinite() || !params.Pos.Rz.allFinite() ||
            !params.Pos.Rr.allFinite() || !params.Pos.Ro.allFinite())
        {
            throw std::invalid_argument("field positions must be finite");
        }
        if (params.Pos.Ro.size() != 0 && params.Pos.Ro.size() != 1 &&
            params.Pos.Ro.size() != params.Pos.Rz.size())
        {
            throw std::invalid_argument(
                "receiver range offsets must contain one value or match receiver depths");
        }
        if (params.MLimit < 1)
        {
            throw std::invalid_argument("mode limit must be positive");
        }
        if (params.SBP.isSet)
        {
            if (params.SBP.theta.size() < 2 ||
                params.SBP.theta.size() != params.SBP.pat.size() ||
                !params.SBP.theta.allFinite() || !params.SBP.pat.allFinite() ||
                (params.SBP.theta.tail(params.SBP.theta.size() - 1) -
                 params.SBP.theta.head(params.SBP.theta.size() - 1)).minCoeff() <= 0.0)
            {
                throw std::invalid_argument(
                    "source beam pattern must be finite and strictly ordered");
            }
        }
    }
}

void validatePublicExecutionCapabilities(const OOKC_parameters &params,
                                         Run_Mode requestedRun)
{
    if (params.freqinfo.Nfreq > 1 || params.freqinfo.freqvec.size() > 1)
    {
        throw std::invalid_argument("Krakenc currently accepts one frequency per solve");
    }
    for (const ssp::Range_Independent_Area &area : params.sspInput)
    {
        static_cast<void>(interpolation(area.SSPType));
        static_cast<void>(boundaryType(area.HSTop.BC));
        static_cast<void>(boundaryType(area.HSBot.BC));
    }
    if (requestedRun != Run_Mode::MODE_M_Modes &&
        params.Pos.GridType != Grid_Mode::MODE_R_Rectangular)
    {
        throw std::invalid_argument("irregular field grid is not implemented");
    }
}

void validatePublicParameters(const OOKC_parameters &params, Run_Mode requestedRun)
{
    validatePublicParameterSemantics(params, requestedRun);
    validatePublicExecutionCapabilities(params, requestedRun);
}

std::vector<AcousticCase> toAcousticCases(const OOKC_parameters &params)
{
    validatePublicParameters(params, Run_Mode::MODE_M_Modes);
    const double frequency = params.freqinfo.freqvec.size() == 1
                                 ? params.freqinfo.freqvec[0]
                                 : params.freqinfo.freq;
    std::vector<AcousticCase> result;
    result.reserve(params.sspInput.size());
    for (const ssp::Range_Independent_Area &area : params.sspInput)
    {
        AcousticCase item;
        item.title = area.Title.empty() ? params.Title : area.Title;
        item.frequency = frequency;
        item.referenceFrequency = params.AttenUnit.referenceFrequency > 0.0
                                      ? params.AttenUnit.referenceFrequency
                                      : frequency;
        item.cLow = params.cLow;
        item.cHigh = params.cHigh;
        item.rMaxKm = params.Rmax;
        item.attenuationUnit = attenuationCode(params.AttenUnit.attnUnit);
        item.absorptionModel = params.AttenUnit.absModel;
        item.volumeAbsorption = params.AttenUnit.volume;
        item.enableRootRestarts = area.enableRootRestarts;
        item.interpolation = interpolation(area.SSPType);
        item.top = makeBoundary(area.HSTop, params.ReflectionCoef.RTop, nullptr,
                                params.ReflectionCoef.isDeg);
        item.bottom = makeBoundary(area.HSBot, params.ReflectionCoef.RBot,
                                   &params.ReflectionCoef.IRC,
                                   params.ReflectionCoef.isDeg);
        const Position &modePosition = params.hasModePos ? params.ModePos : params.Pos;
        item.sourceDepths = toVector(modePosition.Sz);
        item.receiverDepths = toVector(modePosition.Rz);
        for (const ssp::SSPLayer &source : area.layers)
        {
            AcousticLayer layer;
            layer.baseMesh = source.nmesh;
            layer.topDepth = source.z[0];
            layer.bottomDepth = source.z[source.z.size() - 1];
            layer.roughnessRms = source.sigma;
            layer.attenuationPower = source.beta;
            layer.transitionFrequency = source.ft;
            for (Eigen::Index i = 0; i < source.z.size(); ++i)
            {
                layer.samples.push_back({source.z[i], source.alphaR[i], source.betaR[i],
                                         source.rho[i], source.alphaI[i], source.betaI[i]});
            }
            item.layers.push_back(std::move(layer));
        }
        result.push_back(std::move(item));
    }
    return result;
}

FieldParameters toFieldParameters(const OOKC_parameters &params)
{
    if (!params.sspInput.empty())
    {
        validatePublicParameterSemantics(params, Run_Mode::MODE_F_Field);
    }
    validatePublicExecutionCapabilities(params, Run_Mode::MODE_F_Field);
    FieldParameters result;
    result.title = params.Title;
    result.sourceType = sourceCode(params.SourceType);
    result.propagationType = params.modeType == ModeType::Couple ? 'C' : 'A';
    result.coherent = params.coherenceType == CoherenceType::Coherent;
    result.beamPattern = params.SBP.isSet;
    result.modeLimit = params.MLimit;
    result.sourceDepths = toVector(params.Pos.Sz);
    result.receiverDepths = toVector(params.Pos.Rz);
    result.rangesMetres = toVector(params.Pos.Rr);
    if (params.Pos.Ro.size() == params.Pos.Rz.size())
    {
        result.rangeOffsets = toVector(params.Pos.Ro);
    }
    else
    {
        const double offset = params.Pos.Ro.size() == 1 ? params.Pos.Ro[0] : 0.0;
        result.rangeOffsets.assign(result.receiverDepths.size(), offset);
    }
    if (params.RProf.size() != 0)
    {
        for (double value : toVector(params.RProf)) result.profileRangesKm.push_back(value);
    }
    else
    {
        for (const auto &area : params.sspInput) result.profileRangesKm.push_back(area.Range);
    }
    result.beamAnglesDegrees = toVector(params.SBP.theta);
    result.beamAmplitudes = toVector(params.SBP.pat);
    return result;
}

void updatePublicParameters(const std::vector<AcousticCase> &cases,
                            OOKC_parameters &params)
{
    if (cases.empty()) throw std::invalid_argument("acoustic case list must not be empty");
    const AcousticCase &first = cases.front();
    params.Title = first.title;
    params.freqinfo.freq = first.frequency;
    params.freqinfo.Nfreq = 1;
    params.freqinfo.freqvec = Eigen::VectorXd::Constant(1, first.frequency);
    params.cLow = first.cLow;
    params.cHigh = first.cHigh;
    params.Rmax = first.rMaxKm;
    params.AttenUnit.attnUnit = attenuationUnit(first.attenuationUnit);
    params.AttenUnit.absModel = first.absorptionModel;
    params.AttenUnit.referenceFrequency = first.referenceFrequency;
    params.AttenUnit.volume = first.volumeAbsorption;
    params.Pos.Sz = toEigen(first.sourceDepths);
    params.Pos.NSz = static_cast<int>(first.sourceDepths.size());
    params.Pos.Rz = toEigen(first.receiverDepths);
    params.Pos.NRz = static_cast<int>(first.receiverDepths.size());
    params.ModePos = params.Pos;
    params.ModePos.GridType = Grid_Mode::MODE_R_Rectangular;
    params.hasModePos = true;
    params.sspInput.clear();
    params.sspInput.reserve(cases.size());
    for (const AcousticCase &item : cases)
    {
        ssp::Range_Independent_Area area;
        area.Title = item.title;
        area.SSPType = interpolation(item.interpolation);
        area.enableRootRestarts = item.enableRootRestarts;
        area.HSTop = makeHalfSpace(item.top);
        area.HSBot = makeHalfSpace(item.bottom);
        for (const AcousticLayer &source : item.layers)
        {
            ssp::SSPLayer layer;
            layer.npts = static_cast<int>(source.samples.size());
            layer.nmesh = source.baseMesh;
            layer.sigma = source.roughnessRms;
            layer.beta = source.attenuationPower;
            layer.ft = source.transitionFrequency;
            const bool allCsZero =
                !source.samples.empty() &&
                std::all_of(
                    source.samples.begin(), source.samples.end(),
                    [](const AcousticSample &sample) {
                        return sample.cs == 0.0;
                    });
            const bool allCsPositive =
                !source.samples.empty() &&
                std::all_of(
                    source.samples.begin(), source.samples.end(),
                    [](const AcousticSample &sample) {
                        return sample.cs > 0.0;
                    });
            if (!allCsZero && !allCsPositive)
            {
                throw std::invalid_argument(
                    "internal SSP shear sound speeds must be either all zero or all positive");
            }
            layer.Material = allCsPositive
                                 ? Media_Mode::MODE_E_Elastic
                                 : Media_Mode::MODE_A_Acoustic;
            layer.z.resize(layer.npts);
            layer.alphaR.resize(layer.npts);
            layer.alphaI.resize(layer.npts);
            layer.betaR.resize(layer.npts);
            layer.betaI.resize(layer.npts);
            layer.rho.resize(layer.npts);
            for (Eigen::Index i = 0; i < layer.z.size(); ++i)
            {
                const AcousticSample &sample = source.samples[static_cast<std::size_t>(i)];
                layer.z[i] = sample.depth;
                layer.alphaR[i] = sample.cp;
                layer.alphaI[i] = sample.alphaP;
                layer.betaR[i] = sample.cs;
                layer.betaI[i] = sample.alphaS;
                layer.rho[i] = sample.rho;
            }
            area.layers.push_back(std::move(layer));
        }
        params.sspInput.push_back(std::move(area));
    }
    params.NProf = static_cast<int>(cases.size());
    params.RProf = Eigen::VectorXd::Zero(params.NProf);
    params.ReflectionCoef.RTop.resize(1, static_cast<Eigen::Index>(first.top.reflectionSamples.size()));
    params.ReflectionCoef.isDeg = true;
    for (Eigen::Index i = 0; i < params.ReflectionCoef.RTop.size(); ++i)
    {
        const auto &sample = first.top.reflectionSamples[static_cast<std::size_t>(i)];
        params.ReflectionCoef.RTop[i] = {sample.angleDegrees, sample.magnitude, sample.phaseRadians};
    }
    params.ReflectionCoef.RBot.resize(1, static_cast<Eigen::Index>(first.bottom.reflectionSamples.size()));
    for (Eigen::Index i = 0; i < params.ReflectionCoef.RBot.size(); ++i)
    {
        const auto &sample = first.bottom.reflectionSamples[static_cast<std::size_t>(i)];
        params.ReflectionCoef.RBot[i] = {sample.angleDegrees, sample.magnitude, sample.phaseRadians};
    }
    const auto &internal = first.bottom.internalSamples;
    params.ReflectionCoef.IRC.isSet = !internal.empty();
    params.ReflectionCoef.IRC.freq = first.frequency;
    params.ReflectionCoef.IRC.xTab.resize(static_cast<Eigen::Index>(internal.size()));
    params.ReflectionCoef.IRC.fTab.resize(static_cast<Eigen::Index>(internal.size()));
    params.ReflectionCoef.IRC.gTab.resize(static_cast<Eigen::Index>(internal.size()));
    params.ReflectionCoef.IRC.iTab.resize(static_cast<Eigen::Index>(internal.size()));
    for (Eigen::Index i = 0; i < params.ReflectionCoef.IRC.xTab.size(); ++i)
    {
        params.ReflectionCoef.IRC.xTab[i] = internal[static_cast<std::size_t>(i)].eigenvalue;
        params.ReflectionCoef.IRC.fTab[i] = internal[static_cast<std::size_t>(i)].f;
        params.ReflectionCoef.IRC.gTab[i] = internal[static_cast<std::size_t>(i)].g;
        params.ReflectionCoef.IRC.iTab[i] = internal[static_cast<std::size_t>(i)].power10;
    }
}

void updatePublicParameters(const FieldParameters &field,
                            OOKC_parameters &params)
{
    params.Title = field.title;
    params.SourceType = sourceMode(field.sourceType);
    params.modeType = field.propagationType == 'C' ? ModeType::Couple : ModeType::Adiabatic;
    params.coherenceType = field.coherent ? CoherenceType::Coherent : CoherenceType::Incoherent;
    params.MLimit = field.modeLimit;
    params.Pos.Sz = toEigen(field.sourceDepths);
    params.Pos.NSz = static_cast<int>(field.sourceDepths.size());
    params.Pos.Rz = toEigen(field.receiverDepths);
    params.Pos.NRz = static_cast<int>(field.receiverDepths.size());
    params.Pos.Rr = toEigen(field.rangesMetres);
    params.Pos.NRr = static_cast<int>(field.rangesMetres.size());
    params.Pos.Ro = toEigen(field.rangeOffsets);
    params.Pos.NRo = static_cast<int>(field.rangeOffsets.size());
    params.RProf.resize(static_cast<Eigen::Index>(field.profileRangesKm.size()));
    for (Eigen::Index i = 0; i < params.RProf.size(); ++i)
    {
        params.RProf[i] = field.profileRangesKm[static_cast<std::size_t>(i)];
        if (static_cast<std::size_t>(i) < params.sspInput.size())
        {
            params.sspInput[static_cast<std::size_t>(i)].Range = params.RProf[i];
        }
    }
    params.NProf = static_cast<int>(field.profileRangesKm.size());
    params.SBP.isSet = field.beamPattern;
    params.SBP.theta = toEigen(field.beamAnglesDegrees);
    params.SBP.pat = toEigen(field.beamAmplitudes);
    params.SBP.NSBPPts = static_cast<int>(field.beamAnglesDegrees.size());
}
}
