#include "AttenMod.h"

namespace OpenOceanKraken
{
    double parseAttenuation(double freq, double freq0, double alpha, double ft, double beta, double c, const Atten_Mode &AttenUnit)
    {
        double alphaT = 0.0;
        constexpr double dB_TO_NEPERS = 1.0 / 8.6858896;
        constexpr double dB_KHZ_TO_NEPERS = 1.0 / 8685.8896;
        if (is_N(AttenUnit))
        {
            alphaT = alpha;
        }
        else if (is_M(AttenUnit))
        {
            alphaT = alpha * dB_TO_NEPERS;
        }
        else if (is_m(AttenUnit))
        {
            alphaT = alpha * dB_TO_NEPERS;
            if (freq < ft)
                alphaT *= std::pow(freq / freq0, beta);
            else
                alphaT *= (freq / freq0) * std::pow(ft / freq0, beta - 1);
        }
        else if (is_F(AttenUnit))
        {
            alphaT = alpha * freq * dB_KHZ_TO_NEPERS;
        }
        else if (is_W(AttenUnit))
        {
            alphaT = (c != 0.0) ? alpha * freq * dB_TO_NEPERS / c : 0.0;
        }
        else if (is_Q(AttenUnit))
        {
            double omega = 2.0 * pi * freq;
            alphaT = (c * alpha != 0.0) ? omega / (2.0 * c * alpha) : 0.0;
        }
        else if (is_L(AttenUnit))
        {
            double omega = 2.0 * pi * freq;
            alphaT = (c != 0.0) ? alpha * omega / c : 0.0;
        }

        return alphaT;
    }
    double addOceanAbsorption(double alphaT, double freq, const Atten_Mode &AttenUnit)
    {
        if (isThorpe(AttenUnit))
        {
            double f_khz = freq / 1000.0;
            double f2 = f_khz * f_khz;
            double Thorp = (3.3e-3 + 0.11 * f2 / (1.0 + f2) +
                            44.0 * f2 / (4100.0 + f2) + 3e-4 * f2) /
                           8685.8896;
            alphaT += Thorp;
        }
        else if (isFranc_Garr(AttenUnit))
        {
            double FG = Franc_Garr(freq / 1000.0) / 8685.8896;
            alphaT += FG;
        }
        return alphaT;
    }

    std::complex<double> CRCI(double &z, double &c, double &alpha, double &freq, double &freq0,
                              const Atten_Mode &AttenUnit, double &beta, double &ft)
    {

        // 1. 解析用户输入 → 标准衰减
        double alphaT = parseAttenuation(freq, freq0, alpha, ft, beta, c, AttenUnit);

        // 2. 叠加海洋吸收模型
        alphaT = addOceanAbsorption(alphaT, freq, AttenUnit);

        // 3. 转为复声速虚部
        double omega = 2.0 * pi * freq;
        alphaT = (omega != 0.0) ? alphaT * c * c / omega : 0.0;

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
        return AttenUnit.absModel == OceanAbsorptionModel::Thorpe;
    }

    bool isFranc_Garr(const Atten_Mode &AttenUnit) // 判断是否为Franc-Garr沉积层模型
    {
        return AttenUnit.absModel == OceanAbsorptionModel::FrancGarr;
    }

    bool is_F(const Atten_Mode &AttenUnit)
    {
        return AttenUnit.attnUnit == AttenuationUnit::MODE_F_dB_per_m_kHz;
    }
    bool is_L(const Atten_Mode &AttenUnit)
    {
        return AttenUnit.attnUnit == AttenuationUnit::MODE_L_params_lose;
    }
    bool is_M(const Atten_Mode &AttenUnit)
    {
        return AttenUnit.attnUnit == AttenuationUnit::MODE_M_dB_per_m;
    }

    bool is_m(const Atten_Mode &AttenUnit)
    {
        return AttenUnit.attnUnit == AttenuationUnit::MODE_m_dB_per_m;
    }

    bool is_N(const Atten_Mode &AttenUnit)
    {
        return AttenUnit.attnUnit == AttenuationUnit::MODE_N_Nepers_per_m;
    }
    bool is_Q(const Atten_Mode &AttenUnit)
    {
        return AttenUnit.attnUnit == AttenuationUnit::MODE_Q_Quality_Factor;
    }
    bool is_W(const Atten_Mode &AttenUnit)
    {
        return AttenUnit.attnUnit == AttenuationUnit::MODE_W_db_per_lambda;
    }
}
