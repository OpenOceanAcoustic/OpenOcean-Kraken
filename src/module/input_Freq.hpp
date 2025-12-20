#include "paramsBase.h"

class input_Freq : public paramsBase
{
public:
    input_Freq() {}
    virtual ~input_Freq() {}

    virtual void Init(parameters &params) const override 
    {
    }

    virtual void Default(parameters &params) const override
    {
        auto &freqinfo = params.freqinfo;
        freqinfo->Nfreq = 1;
        freqinfo->freq = 200;
        freqinfo->freqvec.resize(1);
        freqinfo->freqvec[0] = freqinfo->freq;
        params.Title = "Pekeris";
        params.RProf.resize(params.NProf);
        params.RProf[0] = 0;
        params.SourceType = Source_Mode::MODE_R_Point;
        params.runMode = Run_Mode::MODE_B_Both;
        params.AttenUnit = Atten_Mode::MODE_W_db_per_lambda;
        params.cLow = 0.0;
        params.cHigh = 20000;
        params.Rmax = 0;
    }

    virtual void Preprocess(parameters &params) const override
    {
        // std::cout<<"input_Freq::Preprocess() is called!"<<std::endl;
    }

    // 添加标题设置函数
    void set_title(parameters &params, const std::string &title)
    {
        auto &titleinfo = params.Title;
        titleinfo = title;
    }

    void set_freq(parameters &params, double freq)
    {
        auto &freqinfo = params.freqinfo;
        freqinfo->Nfreq = 1;
        freqinfo->freq = freq;
    }

    void set_freqvec(parameters &params, VectorXd freqvec)
    {
        auto &freqinfo = params.freqinfo;
        freqinfo->Nfreq = freqvec.size();
        freqinfo->freqvec = freqvec;
    }
    
    // 手动设置全部距离剖面
    void set_RProf(parameters &params, const VectorXd &RProf)
    {
        params.NProf = RProf.size();
        params.RProf = RProf;
    }

    // 插值生成距离剖面
    void set_RProf(parameters &params, const double &start, const double &end, const int &NProf)
    {
        params.RProf.resize(params.NProf);
        linspace(start, end, params.NProf, params.RProf); // 线性插值
    }


    void set_SourceType(parameters &params, Source_Mode type)
    {
        params.SourceType = type;
    }

    void set_RunMode(parameters &params, Run_Mode mode)
    {
        params.runMode = mode;
    }

    void set_cPhase(parameters &params, double cLow, double cHigh)
    {
        params.cLow = cLow;
        params.cHigh = cHigh;
    }

    double get_freq(const parameters &params) const
    {
        auto &freqinfo = params.freqinfo;
        return freqinfo->freq;
    }

    // 打印
    virtual void Echo(parameters &params) const override
    {
    }
};