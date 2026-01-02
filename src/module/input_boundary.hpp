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
        setMunk(params);
    }

    virtual void Preprocess(parameters &params) const override
    {
        double freq = params.freqinfo->freq;
        for (size_t iprof = 0; iprof < params.NProf; iprof++)
        {
            UpdateHSLoss(freq, freq, params.AttenUnit, params.HSTop[iprof], params.HSBot[iprof]);
        }
    }

    // 顶部边界类型选项
    void set_surface_Type(parameters &params, BC_Mode bc, size_t iprof)
    {
        params.HSTop[iprof].BC = bc;
    }
    // 底部边界类型选项
    void set_bottom_Type(parameters &params, BC_Mode bc, size_t iprof)
    {
        params.HSBot[iprof].BC = bc;
    }

    // 设置海底半空间
    void setBottomLine(parameters &params, double zTemp, double alphaR, double alphaI, double betaR, double betaI, double rho, size_t iprof)
    {
        params.HSBot[iprof].Depth = zTemp;
        params.HSBot[iprof].alphaR = alphaR;
        params.HSBot[iprof].alphaI = alphaI;
        params.HSBot[iprof].betaR = betaR;
        params.HSBot[iprof].betaI = betaI;
        params.HSBot[iprof].rho = rho;
    }

    // 设置海面半空间
    void setSurfaceLine(parameters &params, double zTemp, double alphaR, double alphaI, double betaR, double betaI, double rho, size_t iprof)
    {
        params.HSTop[iprof].Depth = zTemp;
        params.HSTop[iprof].alphaR = alphaR;
        params.HSTop[iprof].alphaI = alphaI;
        params.HSTop[iprof].betaR = betaR;
        params.HSTop[iprof].betaI = betaI;
        params.HSTop[iprof].rho = rho;
    }

    // get
    BC_Mode get_surface_Type(parameters &params, size_t iprof)
    {
        return params.HSTop[iprof].BC;
    }

    BC_Mode get_bottom_Type(parameters &params, size_t iprof)
    {
        return params.HSBot[iprof].BC;
    }

    VectorXd get_BottomLine(const parameters &params, size_t iprof)
    {
        VectorXd bottomLine(6);
        bottomLine(0) = params.HSBot[iprof].Depth;
        bottomLine(1) = params.HSBot[iprof].alphaR;
        bottomLine(2) = params.HSBot[iprof].alphaI;
        bottomLine(3) = params.HSBot[iprof].betaR;
        bottomLine(4) = params.HSBot[iprof].betaI;
        bottomLine(5) = params.HSBot[iprof].rho;
        return bottomLine;
    }

    VectorXd get_SurfaceLine(const parameters &params, size_t iprof)
    {
        VectorXd surfaceLine(6);
        surfaceLine(0) = params.HSTop[iprof].Depth;
        surfaceLine(1) = params.HSTop[iprof].alphaR;
        surfaceLine(2) = params.HSTop[iprof].alphaI;
        surfaceLine(3) = params.HSTop[iprof].betaR;
        surfaceLine(4) = params.HSTop[iprof].betaI;
        surfaceLine(5) = params.HSTop[iprof].rho;
        return surfaceLine;
    }

private:
    // 设置默认值的私有方法
    void setDefaultValues(parameters &params) const
    {
        // 第一个剖面 海面、海底参数
        params.HSTop[0].BC = BC_Mode::MODE_V_Vacuum;
        params.HSBot[0].BC = BC_Mode::MODE_A_Half_space;
        params.HSBot[0].alphaI = 0.2;
        params.HSBot[0].alphaR = 1600.0;
        params.HSBot[0].betaI = 0.0;
        params.HSBot[0].betaR = 0.0;
        params.HSBot[0].Depth = 200.0;
        params.HSBot[0].rho = 1.5;
        params.HSBot[0].sigma = 0;
    }

    void setMunk(parameters &params) const
    {
        // 第一个剖面 海面、海底参数
        params.HSTop[0].BC = BC_Mode::MODE_V_Vacuum;
        params.HSBot[0].BC = BC_Mode::MODE_A_Half_space;
        params.HSBot[0].alphaI = 0.2;
        params.HSBot[0].alphaR = 1600.0;
        params.HSBot[0].betaI = 0.0;
        params.HSBot[0].betaR = 0.0;
        params.HSBot[0].Depth = 5000.0;
        params.HSBot[0].rho = 1.5;
        params.HSBot[0].sigma = 0;
    }
    // 打印
    virtual void Echo(parameters &params) const override
    {
    }
};
