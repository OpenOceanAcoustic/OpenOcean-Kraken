#include "module/ParameterAdapters.h"

#include <cassert>
#include <cmath>
#include <limits>
#include <stdexcept>

using namespace OpenOceanKrakenc;

#undef assert
#define assert(condition) do { if (!(condition)) throw std::runtime_error("assertion failed: " #condition); } while (false)

namespace
{
ssp::Range_Independent_Area makeProfile(double rangeKm)
{
    ssp::SSPLayer layer;
    layer.npts = 3;
    layer.nmesh = 20;
    layer.z = (Eigen::VectorXd(3) << 0.0, 50.0, 100.0).finished();
    layer.alphaR = (Eigen::VectorXd(3) << 1500.0, 1490.0, 1510.0).finished();
    layer.alphaI = Eigen::VectorXd::Zero(3);
    layer.betaR = Eigen::VectorXd::Zero(3);
    layer.betaI = Eigen::VectorXd::Zero(3);
    layer.rho = Eigen::VectorXd::Ones(3);

    ssp::Range_Independent_Area area;
    area.Range = rangeKm;
    area.SSPType = SSP_Mode::MODE_C_cLinear;
    area.addLayer(layer);
    area.set_Top_type(BC_Mode::MODE_V_Vacuum);
    area.set_Bottom_Line(100.0, 1800.0, 0.0, 0.0, 0.0, 2.0);
    return area;
}
}

int main()
{
    OOKC_parameters params;
    params.Title = "adapter";
    params.freqinfo.freq = 50.0;
    params.freqinfo.Nfreq = 1;
    params.freqinfo.freqvec = Eigen::VectorXd::Constant(1, 50.0);
    params.cLow = 1400.0;
    params.cHigh = 2000.0;
    params.Rmax = 10.0;
    params.Pos.Sz = Eigen::VectorXd::Constant(1, 25.0);
    params.Pos.NSz = 1;
    params.Pos.Rz = (Eigen::VectorXd(2) << 20.0, 80.0).finished();
    params.Pos.NRz = 2;
    params.Pos.Rr = (Eigen::VectorXd(2) << 1000.0, 2000.0).finished();
    params.Pos.NRr = 2;
    params.Pos.Ro = Eigen::VectorXd::Zero(1);
    params.Pos.NRo = 1;
    params.MLimit = 100;
    params.SourceType = Source_Mode::MODE_S_ScaledCylindrical;
    params.modeType = ModeType::Couple;
    params.coherenceType = CoherenceType::Incoherent;
    params.sspInput = {makeProfile(0.0), makeProfile(5.0)};
    params.NProf = 2;
    params.RProf = (Eigen::VectorXd(2) << 0.0, 5.0).finished();

    validatePublicParameters(params, Run_Mode::MODE_B_Both);
    const auto cases = toAcousticCases(params);
    assert(cases.size() == 2);
    assert(cases[0].title == "adapter");
    assert(cases[0].sourceDepths.size() == 1);
    assert(cases[0].receiverDepths.size() == 2);
    assert(std::abs(cases[1].rMaxKm - 10.0) < 1.0e-12);

    OOKC_parameters attenuationParams = params;
    attenuationParams.AttenUnit.attnUnit = AttenuationUnit::MODE_m_dB_per_m;
    attenuationParams.AttenUnit.absModel = OceanAbsorptionModel::FrancGarr;
    attenuationParams.AttenUnit.referenceFrequency = 40.0;
    attenuationParams.AttenUnit.volume.temperatureCelsius = 10.0;
    attenuationParams.AttenUnit.volume.meanDepthMetres = 500.0;
    attenuationParams.sspInput.front().layers.front().beta = 1.7;
    attenuationParams.sspInput.front().layers.front().ft = 80.0;
    attenuationParams.sspInput.front().layers.front().sigma = 0.02;
    attenuationParams.sspInput.front().HSBot.beta = 1.4;
    attenuationParams.sspInput.front().HSBot.ft = 90.0;
    attenuationParams.sspInput.front().HSBot.sigma = 0.03;
    const auto attenuationCases = toAcousticCases(attenuationParams);
    assert(attenuationCases.front().attenuationUnit == 'm');
    assert(attenuationCases.front().absorptionModel ==
           OceanAbsorptionModel::FrancGarr);
    assert(attenuationCases.front().referenceFrequency == 40.0);
    assert(attenuationCases.front().layers.front().attenuationPower == 1.7);
    assert(attenuationCases.front().layers.front().transitionFrequency == 80.0);
    assert(attenuationCases.front().layers.front().roughnessRms == 0.02);
    assert(attenuationCases.front().bottom.roughnessRms == 0.03);
    OOKC_parameters attenuationRoundTrip;
    updatePublicParameters(attenuationCases, attenuationRoundTrip);
    assert(attenuationRoundTrip.AttenUnit.absModel ==
           OceanAbsorptionModel::FrancGarr);
    assert(attenuationRoundTrip.sspInput.front().layers.front().beta == 1.7);

    const FieldParameters field = toFieldParameters(params);
    assert(field.sourceType == 'S');
    assert(field.propagationType == 'C');
    assert(!field.coherent);
    assert(field.profileRangesKm.size() == 2);
    assert(std::abs(field.profileRangesKm[1] - 5.0) < 1.0e-12);
    assert(field.rangesMetres[1] == 2000.0);

    OOKC_parameters roundTrip;
    updatePublicParameters(cases, roundTrip);
    updatePublicParameters(field, roundTrip);
    assert(roundTrip.NProf == 2);
    assert(roundTrip.sspInput.size() == 2);
    assert(roundTrip.Pos.NSz == 1);
    assert(roundTrip.Pos.NRz == 2);
    assert(roundTrip.Pos.NRr == 2);
    assert(std::abs(roundTrip.RProf[1] - 5.0) < 1.0e-12);
    assert(std::abs(roundTrip.sspInput[1].Range - 5.0) < 1.0e-12);
    assert(roundTrip.SourceType == Source_Mode::MODE_S_ScaledCylindrical);
    assert(roundTrip.modeType == ModeType::Couple);

    {
        OOKC_parameters degrees = params;
        degrees.sspInput.front().HSTop.BC = BC_Mode::MODE_F_File;
        // OpenOcean-Kraken compatibility: false means public degrees;
        // true means the phase values have already been converted to radians.
        degrees.ReflectionCoef.isDeg = false;
        degrees.ReflectionCoef.RTop.resize(1, 2);
        degrees.ReflectionCoef.RTop[0] = {30.0, 0.8, 180.0};
        degrees.ReflectionCoef.RTop[1] = {60.0, 0.7, 90.0};
        const auto converted = toAcousticCases(degrees);
        assert(converted.front().top.reflectionSamples.size() == 2);
        assert(std::abs(converted.front().top.reflectionSamples.front().phaseRadians -
                        std::acos(-1.0)) < 1.0e-12);
    }

    const auto expectRejected = [](const OOKC_parameters &candidate,
                                   Run_Mode mode) {
        bool rejected = false;
        try
        {
            validatePublicParameters(candidate, mode);
        }
        catch (const std::invalid_argument &)
        {
            rejected = true;
        }
        assert(rejected);
    };
    const auto expectFieldConversionRejected = [](const OOKC_parameters &candidate) {
        bool rejected = false;
        try
        {
            static_cast<void>(toFieldParameters(candidate));
        }
        catch (const std::invalid_argument &)
        {
            rejected = true;
        }
        assert(rejected);
    };

    OOKC_parameters pureModField = params;
    pureModField.sspInput.clear();
    pureModField.NProf = 0;
    pureModField.RProf.resize(0);
    {
        OOKC_parameters candidate = pureModField;
        candidate.freqinfo.Nfreq = 2;
        candidate.freqinfo.freqvec =
            (Eigen::VectorXd(2) << 50.0, 100.0).finished();
        expectFieldConversionRejected(candidate);
    }
    {
        OOKC_parameters candidate = pureModField;
        candidate.AttenUnit.absModel = OceanAbsorptionModel::Thorpe;
        static_cast<void>(toFieldParameters(candidate));
    }
    {
        OOKC_parameters candidate = pureModField;
        candidate.Pos.GridType = Grid_Mode::MODE_I_Irregular;
        expectFieldConversionRejected(candidate);
    }

    for (SSP_Mode unsupported : {SSP_Mode::MODE_P_cPCHIP,
                                 SSP_Mode::MODE_S_cCubic,
                                 SSP_Mode::MODE_A_Analytic})
    {
        OOKC_parameters candidate = params;
        candidate.sspInput.front().SSPType = unsupported;
        expectRejected(candidate, Run_Mode::MODE_M_Modes);
    }
    {
        OOKC_parameters candidate = params;
        candidate.sspInput.front().HSTop.BC = BC_Mode::MODE_G_Grain;
        expectRejected(candidate, Run_Mode::MODE_M_Modes);
    }
    for (OceanAbsorptionModel supported : {OceanAbsorptionModel::Thorpe,
                                           OceanAbsorptionModel::FrancGarr,
                                           OceanAbsorptionModel::Biological})
    {
        OOKC_parameters candidate = params;
        candidate.AttenUnit.absModel = supported;
        if (supported == OceanAbsorptionModel::Biological)
        {
            candidate.AttenUnit.volume.biologicalLayers.push_back(
                {0.0, 100.0, 1000.0, 5.0, 2.0});
        }
        validatePublicParameters(candidate, Run_Mode::MODE_M_Modes);
    }
    {
        OOKC_parameters candidate = params;
        candidate.Pos.GridType = Grid_Mode::MODE_I_Irregular;
        expectRejected(candidate, Run_Mode::MODE_F_Field);
    }
    {
        OOKC_parameters candidate = params;
        candidate.Pos.Ro = Eigen::VectorXd::Zero(3);
        candidate.Pos.NRo = 3;
        expectRejected(candidate, Run_Mode::MODE_F_Field);
    }
    {
        OOKC_parameters candidate = params;
        candidate.freqinfo.Nfreq = 2;
        candidate.freqinfo.freqvec =
            (Eigen::VectorXd(2) << 50.0, 100.0).finished();
        expectRejected(candidate, Run_Mode::MODE_M_Modes);
    }
    {
        OOKC_parameters candidate = params;
        candidate.sspInput.front().layers.front().Material =
            Media_Mode::MODE_E_Elastic;
        expectRejected(candidate, Run_Mode::MODE_M_Modes);
    }
    {
        OOKC_parameters candidate = params;
        auto &layer = candidate.sspInput.front().layers.front();
        layer.npts = 1;
        layer.z.conservativeResize(1);
        layer.rho.conservativeResize(1);
        layer.alphaR.conservativeResize(1);
        layer.alphaI.conservativeResize(1);
        layer.betaR.conservativeResize(1);
        layer.betaI.conservativeResize(1);
        expectRejected(candidate, Run_Mode::MODE_M_Modes);
    }
    {
        OOKC_parameters candidate = params;
        candidate.sspInput.front().layers.front().z[1] =
            candidate.sspInput.front().layers.front().z[0];
        expectRejected(candidate, Run_Mode::MODE_M_Modes);
    }
    {
        OOKC_parameters candidate = params;
        candidate.sspInput.front().layers.front().alphaR[1] =
            std::numeric_limits<double>::quiet_NaN();
        expectRejected(candidate, Run_Mode::MODE_M_Modes);
    }
    {
        OOKC_parameters candidate = params;
        candidate.sspInput.front().layers.front().rho[1] = 0.0;
        expectRejected(candidate, Run_Mode::MODE_M_Modes);
    }
    {
        OOKC_parameters candidate = params;
        candidate.sspInput.front().layers.front().alphaI[1] = -1.0;
        expectRejected(candidate, Run_Mode::MODE_M_Modes);
    }
    return 0;
}
