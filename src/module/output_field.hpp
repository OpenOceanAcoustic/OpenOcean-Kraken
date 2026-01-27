#include "paramsBase.h"
namespace OpenOceanKraken
{
    class output_Field : public outputBase
    {
    public:
        output_Field() {};
        virtual ~output_Field() {};

        virtual void Init(OOK_output &output) const override
        {
            output.u_AllSources = nullptr;
            output.v_AllSources = nullptr;
            output.h_AllSources = nullptr;
            // std::cout << "output_Field::Init() is called" << std::endl;
        }

        virtual void Preprocess(OOK_parameters &params, OOK_output &output) const override
        {
            auto &pos = params.Pos;
            size_t N = (size_t)pos.NSz * (size_t)pos.NRz_per_range * (size_t)pos.NRr;
            if (!output.u_AllSources)
            {
                output.u_AllSources = new std::complex<float>[N];
            }
            if (params.is_Velocity && !output.v_AllSources)
            {
                output.v_AllSources = new std::complex<float>[N];
                output.h_AllSources = new std::complex<float>[N];
            }
            std::cout << "OpenOcean-Kraken: output field processed" << std::endl;
        }

        virtual void ClearResults(OOK_parameters &params, OOK_output &output) const override
        {
            auto &pos = params.Pos;
            size_t N = (size_t)pos.NSz * (size_t)pos.NRz_per_range * (size_t)pos.NRr;
            if (output.u_AllSources)
            {
                std::fill(output.u_AllSources, output.u_AllSources + N, std::complex<float>(0.0f, 0.0f));
            }
            if (params.is_Velocity)
            {
                if (output.v_AllSources)
                {
                    std::fill(output.v_AllSources, output.v_AllSources + N, std::complex<float>(0.0f, 0.0f));
                }
                if (output.h_AllSources)
                {
                    std::fill(output.h_AllSources, output.h_AllSources + N, std::complex<float>(0.0f, 0.0f));
                }
            }
        }
        virtual void Finalize(OOK_output &output) const override
        {
            if (output.u_AllSources)
            {
                delete[] output.u_AllSources;
                output.u_AllSources = nullptr;
            }
            if (output.v_AllSources)
            {
                delete[] output.v_AllSources;
                output.v_AllSources = nullptr;
            }
            if (output.h_AllSources)
            {
                delete[] output.h_AllSources;
                output.h_AllSources = nullptr;
            }
        }
    };
}