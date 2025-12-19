#include "paramsBase.h"

class input_sbp : public paramsBase
{
public:
    input_sbp() {}
    virtual ~input_sbp() {}
    virtual void Init(parameters &params) const override
    {
    }

    virtual void Default(parameters &params) const override
    {
        auto &sbp = params.SBP;
        sbp->NSBPPts = 2;
        sbp->theta.resize(sbp->NSBPPts);
        sbp->pat.resize(sbp->NSBPPts);
        sbp->theta(0) = -180.0;
        sbp->theta(1) = 180.0;
        sbp->pat(0) = 1.0;
        sbp->pat(1) = 1.0;
        sbp->isSet = false; // 是否设置了指向性图
    }

    virtual void Preprocess(parameters &params) const override
    {
        auto &sbp = params.SBP;
        if (sbp->isSet)
        {
        }

        // 转化分贝值为标量值
        //  for(int i=0; i<sbp->NSBPPts; i++)
        //  {
        //      sbp->pat(i) = pow(10,sbp->pat(i)/20.0);//convert from dB to linear scale
        //  }
        // std::cout<<"input_sbp::Preprocess() is called!"<<std::endl;
    }

    void set_Pat(parameters &params, const VectorXd &pat, const VectorXd &theta)
    {
        auto &sbp = params.SBP;
        sbp->NSBPPts = pat.size();
        // if(pat.size() == sbp->NSBPPts == theta.size())
        // {
        sbp->theta = theta;
        sbp->pat = pat;
        // 打印 theta 和 pat 的内容
        // std::cout << "theta: " << sbp->theta.transpose() << std::endl;
        // std::cout << "pat: " << sbp->pat.transpose() << std::endl;
        // }
        sbp->isSet = true;
    }

    std::pair<VectorXd, VectorXd> get_SBP(const parameters &params)
    {
        const auto &sbp = params.SBP;
        if (sbp && sbp->NSBPPts > 0 && sbp->theta.size() == sbp->pat.size())
        {
            return {sbp->theta, sbp->pat};
        }
        else
        {
            // 可以返回空向量，或者抛出异常，取决于需求
            return {VectorXd(), VectorXd()};
        }
    }

    // 打印
    virtual void Echo(parameters &params) const override
    {
    }
};
