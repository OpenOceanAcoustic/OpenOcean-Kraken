#include "AttenMod.h"

std::complex<double> CRCI(double &z, double &c, double &alpha, double &freq, double &freq0,
                          const Atten_Mode &AttenUnit, double &beta, double &ft)
{

    double omega = 2.0 * pi * freq;
    double alphaT = 0.0;
    double f2, Thorp, a, FG;
    constexpr float betaPowerLaw = float(1.0);

    if (is_N(AttenUnit))
    {
        alphaT = alpha;
    }
    else if (is_M(AttenUnit))
    {
        alphaT = alpha / 8.6858896;
    }
    else if (is_m(AttenUnit))
    {
        alphaT = alpha / 8.6858896;
        if (freq < ft) // frequency raised to the power beta
            alphaT = alphaT * pow(freq / freq0, beta);
        else // linear in frequency
            alphaT = alphaT * (freq / freq0) * pow(ft / freq0, beta - 1);
    }
    else if (is_F(AttenUnit))
    {
        alphaT = alpha * freq / 8685.8896;
    }
    else if (is_W(AttenUnit))
    {
        if (c != 0.0)
        {
            alphaT = alpha * freq / (8.6858896 * c);
        }
    }
    else if (is_Q(AttenUnit))
    {
        if (c * alpha != 0.0)
        {
            alphaT = omega / (2.0 * c * alpha);
        }
    }
    else if (is_L(AttenUnit))
    {
        if (c != 0.0)
        {
            alphaT = alpha * omega / c;
        }
    }

    if (isThorpe(AttenUnit))
    {
        f2 = pow(freq / 1000.0, 2);
        Thorp = 3.3e-3 + 0.11 * f2 / (1.0 + f2) + 44.0 * f2 / (4100.0 + f2) + 3e-4 * f2;
        Thorp = Thorp / 8685.8896;
        alphaT = alphaT + Thorp;
    }
    else if (isFranc_Garr(AttenUnit))
    {
        // Function Franc_Garr(freq) needs to be defined separately
        FG = Franc_Garr(freq / 1000);
        FG = FG / 8685.8896;
        alphaT = alphaT + FG;
    }

    alphaT = alphaT * c * c / omega;

    // if (alphaT > c) {
    //     PRTFile << "Complex sound speed: " << CRCI << std::endl;
    //     PRTFile << "Usually this means you have an attenuation that is way too high" << std::endl;
    //     // Call ERROUT function with appropriate parameters
    //     ERROUT("AttenMod : CRCI ", "The complex sound speed has an imaginary part > real part");
    // }

    return std::complex<double>(c, alphaT);
}

double Franc_Garr(double f)
{
    double T = 20, Salinity = 35, pH = 8, z_bar = 0;
    const double c = 1412 + 3.21 * T + 1.19 * Salinity + 0.0167 * z_bar;
    double A1, A2, A3, P1, P2, P3, f1, f2;

    // Boric acid contribution
    A1 = 8.86 / c * pow(10, 0.78 * pH - 5);
    P1 = 1;
    f1 = 2.8 * sqrt(Salinity / 35) * pow(10, 4 - 1245 / (T + 273));

    // Magnesium sulfate contribution
    A2 = 21.44 * Salinity / c * (1 + 0.025 * T);
    P2 = 1 - 1.37e-4 * z_bar + 6.2e-9 * pow(z_bar, 2);
    f2 = 8.17 * pow(10, 8 - 1990 / (T + 273)) / (1 + 0.0018 * (Salinity - 35));

    // Viscosity
    P3 = 1 - 3.83e-5 * z_bar + 4.9e-10 * pow(z_bar, 2);
    if (T < 20)
    {
        A3 = 4.937e-4 - 2.59e-5 * T + 9.11e-7 * pow(T, 2) - 1.5e-8 * pow(T, 3);
    }
    else
    {
        A3 = 3.964e-4 - 1.146e-5 * T + 1.45e-7 * pow(T, 2) - 6.5e-10 * pow(T, 3);
    }

    double alpha = A1 * P1 * (f1 * pow(f, 2)) / (pow(f1, 2) + pow(f, 2)) +
                   A2 * P2 * (f2 * pow(f, 2)) / (pow(f2, 2) + pow(f, 2)) +
                   A3 * P3 * pow(f, 2);

    return alpha;
}

bool isThorpe(const Atten_Mode &AttenUnit)
{
    return AttenUnit == Atten_Mode::MODE_FT_dB_per_m_kHz ||
           AttenUnit == Atten_Mode::MODE_LT_params_lose ||
           AttenUnit == Atten_Mode::MODE_MT_dB_per_m ||
           AttenUnit == Atten_Mode::MODE_NT_Nepers_per_m ||
           AttenUnit == Atten_Mode::MODE_QT_Quality_Factor ||
           AttenUnit == Atten_Mode::MODE_WT_db_per_lambda;
}

bool isFranc_Garr(const Atten_Mode &AttenUnit) // 判断是否为Franc-Garr沉积层模型
{
    return AttenUnit == Atten_Mode::MODE_FF_dB_per_m_kHz ||
           AttenUnit == Atten_Mode::MODE_LF_params_lose ||
           AttenUnit == Atten_Mode::MODE_MF_dB_per_m ||
           AttenUnit == Atten_Mode::MODE_NF_Nepers_per_m ||
           AttenUnit == Atten_Mode::MODE_QF_Quality_Factor ||
           AttenUnit == Atten_Mode::MODE_WF_db_per_lambda;
}

bool is_F(const Atten_Mode &AttenUnit)
{
    return AttenUnit == Atten_Mode::MODE_F_dB_per_m_kHz ||
           AttenUnit == Atten_Mode::MODE_FF_dB_per_m_kHz ||
           AttenUnit == Atten_Mode::MODE_FT_dB_per_m_kHz;
}
bool is_L(const Atten_Mode &AttenUnit)
{
    return AttenUnit == Atten_Mode::MODE_L_params_lose ||
           AttenUnit == Atten_Mode::MODE_LF_params_lose ||
           AttenUnit == Atten_Mode::MODE_LT_params_lose;
}
bool is_M(const Atten_Mode &AttenUnit)
{
    return AttenUnit == Atten_Mode::MODE_M_dB_per_m ||
           AttenUnit == Atten_Mode::MODE_MF_dB_per_m ||
           AttenUnit == Atten_Mode::MODE_MT_dB_per_m;
}

bool is_m(const Atten_Mode &AttenUnit)
{
    return AttenUnit == Atten_Mode::MODE_m_dB_per_m ||
           AttenUnit == Atten_Mode::MODE_mF_dB_per_m ||
           AttenUnit == Atten_Mode::MODE_mT_dB_per_m;
}

bool is_N(const Atten_Mode &AttenUnit)
{
    return AttenUnit == Atten_Mode::MODE_N_Nepers_per_m ||
           AttenUnit == Atten_Mode::MODE_NF_Nepers_per_m ||
           AttenUnit == Atten_Mode::MODE_NT_Nepers_per_m;
}
bool is_Q(const Atten_Mode &AttenUnit)
{
    return AttenUnit == Atten_Mode::MODE_Q_Quality_Factor ||
           AttenUnit == Atten_Mode::MODE_QF_Quality_Factor ||
           AttenUnit == Atten_Mode::MODE_QT_Quality_Factor;
}
bool is_W(const Atten_Mode &AttenUnit)
{
    return AttenUnit == Atten_Mode::MODE_W_db_per_lambda ||
           AttenUnit == Atten_Mode::MODE_WF_db_per_lambda ||
           AttenUnit == Atten_Mode::MODE_WT_db_per_lambda;
}