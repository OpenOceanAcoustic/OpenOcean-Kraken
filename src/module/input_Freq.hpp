#include "paramsBase.h"

class input_Freq : public paramsBase
{
public:
    input_Freq() {}
    virtual ~input_Freq() {}

    virtual void Init(parameters &params) const override {}

    virtual void Default(parameters &params) const override
    {
        auto &freqinfo = params.freqinfo;
        freqinfo->Nfreq = 1;
        freqinfo->freq = 100;
        freqinfo->freqvec.resize(1);
        freqinfo->freqvec[0] = 100;

        // std::cout<<"input_Freq::Default() is called!"<<std::endl;
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