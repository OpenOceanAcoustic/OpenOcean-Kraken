
#include "paramsBase.h"

class input_SSP : public paramsBase
{
public:
    input_SSP() {}
    virtual ~input_SSP() {}

    // 初始化函数
    virtual void Init(parameters &params) const override
    {
    }

    // 设置默认值
    virtual void Default(parameters &params) const override
    {
    }

    // 运行前的预处理
    virtual void Preprocess(parameters &params) const override
    {
    }

    // 设置SSP的方法
    void set_SSP(parameters &params, SSP_1D &sspInput, size_t imedia, size_t iprof) const
    {
        auto &ssp = params.SSP[iprof];
        ssp.NPts[imedia] = sspInput.NPts;
        ssp.beta[imedia] = sspInput.beta;
        ssp.ft[imedia] = sspInput.ft;
        ssp.sigma[imedia] = sspInput.sigma;
        // 初始化偏移量
        if (imedia == 0) // 第一个介质层，偏移量为0
        {
            ssp.offset[imedia] = 0;
        }
        else
        {
            ssp.offset[imedia] = ssp.offset[imedia-1] + sspInput.NPts;
        }
        for (int i = 0; i < sspInput.NPts; i++)
        {
            ssp.z[ssp.offset[imedia] + i] = sspInput.z[i];
            ssp.rho[ssp.offset[imedia] + i] = sspInput.rho[i];
            ssp.alphaR[ssp.offset[imedia] + i] = sspInput.alphaR[i];
            ssp.alphaI[ssp.offset[imedia] + i] = sspInput.alphaI[i];
            ssp.betaR[ssp.offset[imedia] + i] = sspInput.betaR[i];
            ssp.betaI[ssp.offset[imedia] + i] = sspInput.betaI[i];
        }
    }

    void set_AttenUnit(parameters &params, Atten_Mode unit) const
    {
    }

    // 设置SSP模式
    void set_SSP_mode(parameters &params, SSP_Mode mode) const
    {
    }

    SSP_1D get_SSP(const parameters &params, size_t imedia, size_t iprof) const
    {
        auto &ssp = params.SSP[iprof];
        SSP_1D sspOutput;
        int NPts = ssp.NPts[imedia];
        int offset = ssp.offset[imedia];
        sspOutput.NPts = NPts;
        sspOutput.beta = ssp.beta[imedia];
        sspOutput.ft = ssp.ft[imedia];
        sspOutput.sigma = ssp.sigma[imedia];
        sspOutput.z = ssp.z.segment(offset, NPts);
        sspOutput.rho = ssp.rho.segment(offset, NPts);
        sspOutput.alphaR = ssp.alphaR.segment(offset, NPts);
        sspOutput.alphaI = ssp.alphaI.segment(offset, NPts);
        sspOutput.betaR = ssp.betaR.segment(offset, NPts);
        sspOutput.betaI = ssp.betaI.segment(offset, NPts);
    }

    // 打印
    virtual void Echo(parameters &params) const override
    {
    }

private:
    // 设置默认值的私有方法
    void setDefaultValues(parameters &params) const
    {
    }
};
