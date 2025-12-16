#include "paramsBase.h"
// 描述地形的
class input_Boundary : public paramsBase
{
public:
    input_Boundary() {}
    virtual ~input_Boundary() {}

    virtual void Init(parameters &params) const override
    {
    }

    virtual void Default(parameters &params) const override
    {
    }

    virtual void Preprocess(parameters &params) const override
    {
        if (!params.isBotSet)
            Default_Bot(params);
        if (!params.isTopSet)
            Default_Top(params);

        // constexpr double betaPowerLaw = 1.0;
        ComputeBdryTangentNormal(params);
        // 对BC也要进行处理
        auto &HS = params.Bdry->Top.HS;
        switch (HS.BC)
        {
        case BC_Mode::MODE_A_Half_space:
            // 计算cP和cS
            HS.cP = CRCI(HS.Depth,
                         HS.alphaR,
                         HS.alphaI,
                         params.freqinfo->freq,
                         params.freqinfo->freq,
                         params.SSP->AttenUnit,
                         params.ft);

            HS.cS = CRCI(HS.Depth,
                         HS.betaR,
                         HS.betaI,
                         params.freqinfo->freq,
                         params.freqinfo->freq,
                         params.SSP->AttenUnit,
                         params.ft);

            break;

        case BC_Mode::MODE_G_Grain:
            double vr, rhor, alpha2_f;
            double &Mz = HS.Mz;
            if (Mz >= -1 && Mz < 1)
            {
                vr = 0.002709 * pow(Mz, 2) - 0.056452 * Mz + 1.2778;
                rhor = 0.007797 * pow(Mz, 2) - 0.17057 * Mz + 2.3139;
            }
            else if (Mz >= 1 && Mz < 5.3)
            {
                vr = -0.0014881 * pow(Mz, 3) + 0.0213937 * pow(Mz, 2) - 0.1382798 * Mz + 1.3425;
                rhor = -0.0165406 * pow(Mz, 3) + 0.2290201 * pow(Mz, 2) - 1.1069031 * Mz + 3.0455;
            }
            else
            {
                vr = -0.0024324 * Mz + 1.0019;
                rhor = -0.0012973 * Mz + 1.1565;
            }

            if (Mz >= -1 && Mz < 0)
            {
                alpha2_f = 0.4556;
            }
            else if (Mz >= 0 && Mz < 2.6)
            {
                alpha2_f = 0.4556 + 0.0245 * Mz;
            }
            else if (Mz >= 2.6 && Mz < 4.5)
            {
                alpha2_f = 0.1978 + 0.1245 * Mz;
            }
            else if (Mz >= 4.5 && Mz < 6.0)
            {
                alpha2_f = 8.0399 - 2.5228 * Mz + 0.20098 * pow(Mz, 2);
            }
            else if (Mz >= 6.0 && Mz < 9.5)
            {
                alpha2_f = 0.9431 - 0.2041 * Mz + 0.0117 * pow(Mz, 2);
            }
            else
            {
                alpha2_f = 0.0601;
            }
            double alphaR = vr * 1500.0;
            double alphaI = alpha2_f * (vr / 1000) * 1500.0 * log(10.0) / (40.0 * pi); // loss parameter Sect. IV., Eq. (4) of handbook
            HS.cP = CRCI(HS.Depth,
                         alphaR,
                         alphaI,
                         params.freqinfo->freq,
                         params.freqinfo->freq,
                         params.BG_AttenUnit,
                         params.ft);

            HS.cS = 0.0;
            HS.rho = rhor;
            break;
        }

        auto &HS1 = params.Bdry->Bot.HS;
        switch (HS1.BC)
        {
        case BC_Mode::MODE_A_Half_space:
            // 计算cP和cS
            HS1.cP = CRCI(HS1.Depth,
                          HS1.alphaR,
                          HS1.alphaI,
                          params.freqinfo->freq,
                          params.freqinfo->freq,
                          params.SSP->AttenUnit,
                          params.ft);

            HS1.cS = CRCI(HS1.Depth,
                          HS1.betaR,
                          HS1.betaI,
                          params.freqinfo->freq,
                          params.freqinfo->freq,
                          params.SSP->AttenUnit,
                          params.ft);

            break;

        case BC_Mode::MODE_G_Grain:
            double vr, rhor, alpha2_f;
            double &Mz = HS1.Mz;
            if (Mz >= -1 && Mz < 1)
            {
                vr = 0.002709 * pow(Mz, 2) - 0.056452 * Mz + 1.2778;
                rhor = 0.007797 * pow(Mz, 2) - 0.17057 * Mz + 2.3139;
            }
            else if (Mz >= 1 && Mz < 5.3)
            {
                vr = -0.0014881 * pow(Mz, 3) + 0.0213937 * pow(Mz, 2) - 0.1382798 * Mz + 1.3425;
                rhor = -0.0165406 * pow(Mz, 3) + 0.2290201 * pow(Mz, 2) - 1.1069031 * Mz + 3.0455;
            }
            else
            {
                vr = -0.0024324 * Mz + 1.0019;
                rhor = -0.0012973 * Mz + 1.1565;
            }

            if (Mz >= -1 && Mz < 0)
            {
                alpha2_f = 0.4556;
            }
            else if (Mz >= 0 && Mz < 2.6)
            {
                alpha2_f = 0.4556 + 0.0245 * Mz;
            }
            else if (Mz >= 2.6 && Mz < 4.5)
            {
                alpha2_f = 0.1978 + 0.1245 * Mz;
            }
            else if (Mz >= 4.5 && Mz < 6.0)
            {
                alpha2_f = 8.0399 - 2.5228 * Mz + 0.20098 * pow(Mz, 2);
            }
            else if (Mz >= 6.0 && Mz < 9.5)
            {
                alpha2_f = 0.9431 - 0.2041 * Mz + 0.0117 * pow(Mz, 2);
            }
            else
            {
                alpha2_f = 0.0601;
            }
            double alphaR = vr * 1500.0;
            double alphaI = alpha2_f * (vr / 1000) * 1500.0 * log(10.0) / (40.0 * pi); // loss parameter Sect. IV., Eq. (4) of handbook
            double ft = 1000;
            HS1.cP = CRCI(HS.Depth,
                          alphaR,
                          alphaI,
                          params.freqinfo->freq,
                          params.freqinfo->freq,
                          params.BG_AttenUnit,
                          ft);

            HS1.cS = 0.0;
            HS1.rho = rhor;
            break;
        }

        // std::cout<<"input_Boundary::Preprocess() is called"<<std::endl;
    }

    // 顶部边界类型选项
    void set_surface_Type(parameters &params, BC_Mode bc)
    {
    }
    // 底部边界类型选项
    void set_bottom_Type(parameters &params, BC_Mode bc)
    {
    }

    // 设置海底半空间
    void setBottomLine(parameters &params, double zTemp, double alphaR, double alphaI, double betaR, double betaI, double rho)
    {
    }
    // 设置海面半空间
    void setSurfaceLine(parameters &params, double zTemp, double alphaR, double alphaI, double betaR, double betaI, double rho)
    {
    }
    // 设置粒子边界条件
    void setSurface_Grain(parameters &params, double zTemp, double Mz)
    {
    }
    // 设置粒子边界条件
    void setBottom_Grain(parameters &params, double zTemp, double Mz)
    {
    }

    // get
    BC_Mode get_surface_Type(parameters &params)
    {
    }

    BC_Mode get_bottom_Type(parameters &params)
    {
    }

private:
    // 打印
    virtual void Echo(parameters &params) const override
    {
    }
};
