#include "paramsBase.h"

class output_Eigen : public outputBase
{
public:
    output_Eigen() {};
    virtual ~output_Eigen() {};

    virtual void Init(kkc_output &output) const override
    {
        output.eigen = nullptr;
    }

    virtual void Preprocess(parameters &params, kkc_output &output) const override
    {
        auto &pos = params.Pos;
        double freq = params.freqinfo->freq;
        if (!output.eigen)
        {
            output.eigen = new EigenParams[params.NProf];
        }
        for (int iprof = 0; iprof < params.NProf; iprof++)
        {
            double cmin = params.SSP[iprof].alphaR.minCoeff();//找元素中最小值
            output.eigen[iprof].firstM = (size_t)(2.0 * params.SSP[iprof].depth.tail(1)(0) * freq / cmin * 1.1 + 10);
            output.eigen[iprof].resize(output.eigen[iprof].firstM, pos->NSz, pos->NRz, params.NMeshMax, params.mesh.NSets);
        }
    }

    virtual void ClearResults(parameters &params, kkc_output &output) const override
    {
        for (int iprof = 0; iprof < params.NProf; iprof++)
        {
            output.eigen[iprof].setZero();
        }
    }
    virtual void Finalize(kkc_output &output) const override
    {
        delete[] output.eigen;
        output.eigen = nullptr;
    }
};
