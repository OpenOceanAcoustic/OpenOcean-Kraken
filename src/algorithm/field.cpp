#include "field.h"

void field(EigenFunction &eigenfun, EigenParams &eigen, parameters &params, std::complex<float> *uAllSources, int isz)
{
    VectorXcd constt, sumk;
    constt.resize(eigen.M);
    sumk = VectorXcd::Zero(eigen.M);

    for (int i = 0; i < eigen.M; i++)
    {
        constt(i) = I1D * std::sqrt(2.0 * pi) * std::exp(I1D * pi / 4.0) * eigenfun.phiS(i, isz);
    }

    VectorXcd Hank(eigen.M);
    // VectorXcd phi(eigen.M);
    double rLeft;

    for (int irr = 0; irr < params.Pos->NRr; irr++)
    {   
        if (irr > 0)
            rLeft = std::max(params.Pos->Rr(irr - 1), 0.0);
        else
            rLeft = 0.0;

        double Rr = params.Pos->Rr(irr);
        for (int i = 0; i < eigen.M; i++)
        {
            sumk(i) = sumk(i) + eigen.k(i) * (Rr - rLeft);

            Hank(i) = constt(i) * std::exp(-I1D * sumk(i));  // coherent   case

            if (params.SourceType == Source_Mode::MODE_R_Point)  // Cylindrical coordinates
            {
                if (Rr == 0.0)
                    Hank(i) = 0.0;
                else
                    Hank(i) = Hank(i) / std::sqrt(eigen.k(i) * Rr);
            }
            else if (params.SourceType == Source_Mode::MODE_X_Line)  // Cartesian coordinates
            {
                Hank(i) = Hank(i) / eigen.k(i);
            }
            else  // Scaled cylindrical coordinates
            {
                Hank(i) = Hank(i) / std::sqrt(eigen.k(i));
            }
        }
        
        // For each receiver, add up modal contributions
        for (int irz = 0; irz < params.Pos->NRz; irz++)
        {
            size_t base = GetFieldAddr(isz, irz, irr, &params.Pos[0]);
            for (int i = 0; i < eigen.M; i++)
            {
                uAllSources[base] += eigenfun.phiR(i, irz) * Hank(i);
            }
        }
    }
}

inline size_t GetFieldAddr(int32_t isz, int32_t id, int32_t ir, const Position *Pos)
{
    return ((size_t)isz * (size_t)Pos->NRz_per_range + (size_t)id) * (size_t)Pos->NRr + (size_t)ir;
}
