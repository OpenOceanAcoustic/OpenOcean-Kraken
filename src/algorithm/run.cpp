#include "run.h"

// 计算本征值和本征函数
void EigenVWorker(size_t iprof, parameters &params, TridMtx& trid, kkc_output &output)
{
    double freq = params.freqinfo->freq;
    auto &ssp = params.SSP[iprof];
    auto HSTop = params.HSTop[iprof];
    auto HSBot = params.HSBot[iprof];
    auto eigen = output.eigen[iprof];
    double error;

    for (int iset = 0; iset < params.mesh.NSets; iset++)
    {
        int ntimes = params.mesh.NV[iset];
        TridPreprocess(iset, iprof, params, trid, ntimes);
        SolveEp(iset, iprof, params.mesh.NSets, eigen, trid, params, error);
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

    size_t m = 0;
    while (m < eigen.M && eigen.Extrap(m) > SQ(2 * pi * params.freqinfo->freq / params.Chigh))
    {
        m++;
    }
    eigen.M = m;

    for (int i = 0; i < m; i++)
    {
        eigen.k(i) = sqrt(eigen.Extrap(i) + eigen.k(i));
    }
}

// 计算声压
void ComputePressure(const size_t& iprof, parameters &params, kkc_output &output)
{
    auto eigen = output.eigen[iprof];
    for (int isz = 0; isz < params.Pos->NSz; isz++)
    {
        Evaluate(eigen, params, isz, output.u_AllSources, output.v_AllSources, output.h_AllSources);
    }
}