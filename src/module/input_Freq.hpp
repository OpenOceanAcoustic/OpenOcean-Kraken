#include "paramsBase.h"

namespace OpenOceanKraken
{

    class input_Freq : public paramsBase
    {
    public:
        input_Freq() {}
        virtual ~input_Freq() {}

        virtual void Init(OOK_parameters &params) const override
        {
        }

        virtual void Default(OOK_parameters &params) const override
        {
            setMunk(params);
        }

        virtual void Preprocess(OOK_parameters &params) const override
        {
            // std::cout<<"input_Freq::Preprocess() is called!"<<std::endl;
        }

        // 添加标题设置函数
        void set_title(OOK_parameters &params, const std::string &title)
        {
            auto &titleinfo = params.Title;
            titleinfo = title;
        }

        void set_freq(OOK_parameters &params, double freq)
        {
            auto &freqinfo = params.freqinfo;
            freqinfo.Nfreq = 1;
            freqinfo.freq = freq;
        }

        void set_freqvec(OOK_parameters &params, Eigen::VectorXd freqvec)
        {
            auto &freqinfo = params.freqinfo;
            freqinfo.Nfreq = freqvec.size();
            freqinfo.freqvec = freqvec;
        }


        void set_SourceType(OOK_parameters &params, Source_Mode type)
        {
            params.SourceType = type;
        }

        void set_RunMode(OOK_parameters &params, Run_Mode mode)
        {
            params.runMode = mode;
        }

        void set_cPhase(OOK_parameters &params, double cLow, double cHigh)
        {
            params.cLow = cLow;
            params.cHigh = cHigh;
        }

        void setMunk(OOK_parameters &params) const
        {
            auto &freqinfo = params.freqinfo;
            freqinfo.Nfreq = 1;
            freqinfo.freq = 50;
            freqinfo.freqvec.resize(1);
            freqinfo.freqvec[0] = freqinfo.freq;
            params.Title = "Munk";
            // params.RProf.resize(1);
            // params.RProf[0] = 0;
            params.SourceType = Source_Mode::MODE_R_Point;
            params.runMode = Run_Mode::MODE_B_Both;
            params.AttenUnit = Atten_Mode{};
            params.cLow = 1500;
            params.cHigh = 1600;
            params.Rmax = 50e3;
        }
        void set_Velocity_enable(OOK_parameters &params, bool is_Velocity)
        {
            params.is_Velocity = is_Velocity;
        }


        // 打印
        virtual void Echo(OOK_parameters &params) const override
        {
        }

        

    };
}
