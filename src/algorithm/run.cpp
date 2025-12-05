#include "run.h"

void run()
{
    parameters params;
    set_pekeris(params);
    // set_Munk(params);
    // 初始化cp cs
    UpdateSSPLoss(params.freqinfo->freq, params.freqinfo->freq, params.NMedia,
                  params.SSPType, params.AttenUnit, params.SSP);
    UpdateHSLoss(params.freqinfo->freq, params.freqinfo->freq, params.NMedia,
                 params.AttenUnit, params.HSTop, params.HSBot);

    // TODO 计算本征值和本征函数
    KrakenMatrix kramtrx;
    EigenParams eigen;
    EigenFunction eigenfun;
    double error;
    params.mesh.hV.resize(params.mesh.NSets);
    for (int iset = 0; iset < params.mesh.NSets; iset++)
    {
        int ntimes = params.mesh.NV[iset];
        Initialize(iset, params, kramtrx, ntimes);
        SolveEp(iset, params.mesh.NSets, eigen, eigenfun, kramtrx, params, error);
        if (error * params.Rmax < 1.0)
        {
            break;
        }
        else
        {
            if (iset == params.mesh.NSets - 1)
                cout << "Warning in KRAKEN : Too many meshes needed: check convergence" << endl;
        }
    }

    int M = 0;
    VectorXd Ex1 = eigen.Extrap.segment(0, eigen.M);
    // cout << "Ex1: \n"
    //      << Ex1 << endl;
    while (M < eigen.M && Ex1(M) > SQ(2 * pi * params.freqinfo->freq / params.Chigh))
    {
        M++;
    }
    eigen.M = M;
    eigen.k.resize(M);
    for (int i = 0; i < M; i++)
    {
        eigen.k(i) = sqrt(eigen.Extrap(i) + eigen.k(i));
    }
    size_t N = (size_t)params.Pos->NSz * (size_t)params.Pos->NRz_per_range * (size_t)params.Pos->NRr;
    std::complex<float> *u_AllSources;
    u_AllSources = new std::complex<float>[N];
    std::complex<float> *uAllSources_vr;
    uAllSources_vr = new std::complex<float>[N];
    std::complex<float> *uAllSources_vz;
    uAllSources_vz = new std::complex<float>[N];
    for (int isz = 0; isz < params.Pos->NSz; isz++)
    {
        Evaluate(eigenfun, eigen, params, isz, u_AllSources, uAllSources_vr, uAllSources_vz);
    }
    string filename = "test_pressure";
    string filename_vr = "test_vr";
    string filename_vz = "test_vz";
    export_shd(filename, params, u_AllSources);
    export_shd(filename_vr, params, uAllSources_vr);
    export_shd(filename_vz, params, uAllSources_vz);
}

