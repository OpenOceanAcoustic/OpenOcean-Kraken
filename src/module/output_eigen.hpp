#include "paramsBase.h"
namespace OpenOceanKraken
{

    class output_Eigen : public outputBase
    {
    public:
        output_Eigen() {};
        virtual ~output_Eigen() {};

        virtual void Init(OOK_output &output) const override
        {
            output.eigen = nullptr;
        }

        virtual void Preprocess(OOK_parameters &params, OOK_output &output) const override
        {
            auto &pos = params.Pos;
            double freq = params.freqinfo.freq;
            if (!output.eigen){
                
                output.eigen = new EigenParams[params.SSP.size()];
            }
            for (int iprof = 0; iprof < params.SSP.size(); iprof++)
            {
                double cmin = params.SSP.at(iprof).alphaR.minCoeff();
                output.eigen[iprof].firstM = (size_t)(2.0 * params.SSP.at(iprof).depth.tail(1)(0) * freq / cmin * 1.1 + 10);
                output.eigen[iprof].resize(output.eigen[iprof].firstM, pos.NSz, pos.NRz, params.NMeshMax, params.mesh.NSets);
            }

            std::cout << "OpenOcean-Kraken: output eigen processed" << std::endl;
        }

        virtual void ClearResults(OOK_parameters &params, OOK_output &output) const override
        {
            for (int iprof = 0; iprof < params.SSP.size(); iprof++)
            {
                output.eigen[iprof].setZero();
            }
        }
        virtual void Finalize(OOK_output &output) const override
        {
            delete[] output.eigen;
            output.eigen = nullptr;
        }
    };
}