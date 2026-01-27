#include "paramsBase.h"

namespace OpenOceanKraken
{
    class input_reflcoef : public paramsBase
    {
    public:
        input_reflcoef() {}
        virtual ~input_reflcoef() {}

        virtual void Init(OOK_parameters &params) const override {}

        virtual void Default(OOK_parameters &params) const override
        {
            setMunk(params);

        }

        virtual void Preprocess(OOK_parameters &params) const override
        {
            auto &RTop = params.ReflectionCoef.RTop;
            auto &RBot = params.ReflectionCoef.RBot;
            bool &isDeg = params.ReflectionCoef.isDeg;
            if (isDeg)
            {
                for (int i = 0; i < RTop.size(); ++i)
                {
                    RTop(i).phi *= DegRad;
                }
                for (int i = 0; i < RBot.size(); ++i)
                {
                    RBot(i).phi *= DegRad;
                }
                isDeg = true;
            }

            // std::cout<<"input_reflcoef::Preprocess() is called!"<<std::endl;
        }
        // 设置顶部反射系数
        void set_ReflCoef_Top(OOK_parameters &params, std::vector<ReflectionCoef> ReflCoef)
        {
            auto &RTop = params.ReflectionCoef.RTop;
            RTop.resize(ReflCoef.size());
            for (int i = 0; i < ReflCoef.size(); ++i)
            {
                RTop(i).R = ReflCoef[i].R;
                RTop(i).theta = ReflCoef[i].theta;
                RTop(i).phi = ReflCoef[i].phi;
            }
            params.Bdry.Top.BC = BC_Mode::MODE_F_File; // 边界条件设置为文件读取
        }

        // 设置底部反射系数
        void set_ReflCoef_Bot(OOK_parameters &params, std::vector<ReflectionCoef> ReflCoef)
        {
            auto &RBot = params.ReflectionCoef.RBot;
            RBot.resize(ReflCoef.size());
            for (int i = 0; i < ReflCoef.size(); ++i)
            {
                RBot(i).R = ReflCoef[i].R;
                RBot(i).theta = ReflCoef[i].theta;
                RBot(i).phi = ReflCoef[i].phi;
            }
            params.Bdry.Bot.BC = BC_Mode::MODE_F_File; // 边界条件设置为文件读取
        }

        std::vector<ReflectionCoef> get_ReflCoef_Top(const OOK_parameters &params)
        {
            std::vector<ReflectionCoef> ReflCoef(params.ReflectionCoef.RTop.size());
            bool isDeg = params.ReflectionCoef.isDeg;
            for (int i = 0; i < params.ReflectionCoef.RTop.size(); ++i)
            {
                ReflCoef[i].R = params.ReflectionCoef.RTop(i).R;
                ReflCoef[i].phi = params.ReflectionCoef.RTop(i).phi;
                ReflCoef[i].theta = params.ReflectionCoef.RTop(i).theta;
                if (isDeg)
                {
                    ReflCoef[i].phi *= RadDeg;
                }
            }
            return ReflCoef;
        }

        std::vector<ReflectionCoef> get_ReflCoef_Bot(const OOK_parameters &params)
        {
            std::vector<ReflectionCoef> ReflCoef(params.ReflectionCoef.RBot.size());
            bool isDeg = params.ReflectionCoef.isDeg;
            for (int i = 0; i < params.ReflectionCoef.RBot.size(); ++i)
            {
                ReflCoef[i].R = params.ReflectionCoef.RBot(i).R;
                ReflCoef[i].phi = params.ReflectionCoef.RBot(i).phi;
                ReflCoef[i].theta = params.ReflectionCoef.RBot(i).theta;
                if (isDeg)
                {
                    ReflCoef[i].phi *= RadDeg;
                }
            }
            return ReflCoef;
        }

        void setMunk(OOK_parameters &params) const
        {
            auto &RTop = params.ReflectionCoef.RTop;
            auto &RBot = params.ReflectionCoef.RBot;
            RTop.resize(1);
            RBot.resize(1);
            RTop(0).R = 1.0;
            RTop(0).theta = 0.0;
            RTop(0).phi = 0.0;
            RBot(0).R = 1.0;
            RBot(0).theta = 0.0;
            RBot(0).phi = 0.0;
        }


        // 打印
        virtual void Echo(OOK_parameters &params) const override
        {
        }
    };
}