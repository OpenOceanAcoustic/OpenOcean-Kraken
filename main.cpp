#include "kkc_params.h"
#include "sspMod.h"

int main()
{
    double freq = 1000.0;
    SSPStructure ssp;
    ssp.Type = SSP_Mode::MODE_N_n2Linear;

    ssp.NMedia = 2;
    ssp.Nz = 4;
    ssp.N.resize(ssp.NMedia);
    ssp.L.resize(ssp.NMedia);
    ssp.NPts.resize(ssp.NMedia);
    ssp.Loc.resize(ssp.NMedia);
    ssp.z.resize(ssp.Nz);
    ssp.alphaR.resize(ssp.Nz);
    ssp.betaR.resize(ssp.Nz);
    ssp.rho.resize(ssp.Nz);
    ssp.alphaI.resize(ssp.Nz);
    ssp.betaI.resize(ssp.Nz);
    ssp.beta.resize(ssp.Nz);
    ssp.ft.resize(ssp.Nz);

    ssp.NPts << 2, 2;
    ssp.Loc << 0, 2;
    ssp.N << 10, 10;
    ssp.L << 0, 10;
    ssp.z << 0, 1000.0, 1000.0, 2000;
    ssp.alphaR << 1550, 1530, 1600, 1650;
    ssp.betaR << 0.00, 0.00, 10, 15;
    ssp.rho << 1.03, 1.03, 1.3, 1.6;
    ssp.alphaI << 0.0, 0.0, 0.2, 0.15;
    ssp.betaI << 0.0, 0.0, 0.02, 0.02;
    ssp.beta << 0.001, 0.001, 0.001, 0.001;
    ssp.ft << 1000.0, 1000.0, 1000.0, 1000.0;
    ssp.AttenUnit = Atten_Mode::MODE_W_db_per_lambda;
    ssp.cp.resize(ssp.Nz);
    ssp.cs.resize(ssp.Nz);


    for (size_t iz = 0; iz < ssp.Nz; iz++)
    {
        ssp.cp(iz) = CRCI(ssp.z(iz), ssp.alphaR(iz), ssp.alphaI(iz), freq, freq,
                          ssp.AttenUnit, ssp.beta(iz), ssp.ft(iz));
        ssp.cs(iz) = CRCI(ssp.z(iz), ssp.betaR(iz), ssp.betaI(iz), freq, freq,
                          ssp.AttenUnit, ssp.beta(iz), ssp.ft(iz));
    }

    ssp.cp_int.resize(ssp.N.sum());
    ssp.cs_int.resize(ssp.N.sum());
    ssp.rho_int.resize(ssp.N.sum());

    // 调用EvaluateSSP函数
    EvaluateSSP(ssp);

    // 打印
    std::cout << "cp: " << ssp.cp_int.transpose() << std::endl;
    std::cout << "cs: " << ssp.cs_int.transpose() << std::endl;
    std::cout << "rho_k: " << ssp.rho_int.transpose() << std::endl;
}