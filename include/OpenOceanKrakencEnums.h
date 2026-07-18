#ifndef OPEN_OCEAN_KRAKENC_ENUMS_H
#define OPEN_OCEAN_KRAKENC_ENUMS_H

#include <vector>

namespace OpenOceanKrakenc
{
enum class SSP_Mode
{
    MODE_N_n2Linear,
    MODE_C_cLinear,
    MODE_P_cPCHIP,
    MODE_S_cCubic,
    MODE_A_Analytic
};

enum class Media_Mode
{
    MODE_A_Acoustic,
    MODE_E_Elastic
};

enum class AttenuationUnit
{
    MODE_F_dB_per_m_kHz,
    MODE_L_params_lose,
    MODE_M_dB_per_m,
    MODE_m_dB_per_m,
    MODE_N_Nepers_per_m,
    MODE_Q_Quality_Factor,
    MODE_W_db_per_lambda
};

enum class OceanAbsorptionModel
{
    None,
    Thorpe,
    FrancGarr,
    Biological
};

struct BiologicalAbsorptionLayer
{
    double topDepthMetres = 0.0;
    double bottomDepthMetres = 0.0;
    double resonanceFrequencyHz = 0.0;
    double qualityFactor = 0.0;
    double peakAttenuationDbPerKm = 0.0;
};

struct VolumeAbsorptionParameters
{
    double temperatureCelsius = 20.0;
    double salinityPsu = 35.0;
    double ph = 8.0;
    double meanDepthMetres = 0.0;
    std::vector<BiologicalAbsorptionLayer> biologicalLayers;
};

struct AttenuationContext
{
    char unit = 'W';
    OceanAbsorptionModel model = OceanAbsorptionModel::None;
    double frequency = 0.0;
    double referenceFrequency = 0.0;
    double beta = 0.0;
    double transitionFrequency = 0.0;
    VolumeAbsorptionParameters volume;
};

enum class BC_Mode
{
    MODE_R_Rigid,
    MODE_V_Vacuum,
    MODE_F_File,
    MODE_A_Half_space,
    MODE_G_Grain,
    MODE_P_Precomputed
};

enum class Source_Mode
{
    MODE_R_Point,
    MODE_X_Line,
    MODE_S_ScaledCylindrical
};

enum class Grid_Mode
{
    MODE_R_Rectangular,
    MODE_I_Irregular
};

enum class Run_Mode
{
    MODE_M_Modes,
    MODE_F_Field,
    MODE_B_Both
};

enum class CoherenceType
{
    Coherent,
    Incoherent
};

enum class ModeType
{
    Adiabatic,
    Couple
};

struct Atten_Mode
{
    AttenuationUnit attnUnit = AttenuationUnit::MODE_W_db_per_lambda;
    OceanAbsorptionModel absModel = OceanAbsorptionModel::None;
    double referenceFrequency = 0.0;
    VolumeAbsorptionParameters volume;
};
}

#endif
