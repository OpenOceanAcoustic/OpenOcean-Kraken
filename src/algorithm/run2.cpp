#include "run2.h"

// 计算本征值和本征函数
void EigenVWorker(size_t iprof, parameters &params, kkc_output &output)
{
    double freq = params.freqinfo->freq;
    auto &ssp = params.SSP[iprof];
    auto HSTop = params.HSTop[iprof];
    auto HSBot = params.HSBot[iprof];
    auto eigen = output.eigen[iprof];
    KrakenMatrix kramtrx;
    double error;

    for (int iset = 0; iset < params.mesh.NSets; iset++)
    {
        int ntimes = params.mesh.NV[iset];
        Initialize(iset, iprof, params, kramtrx, ntimes);
        SolveEpMode(iset, iprof, params.mesh.NSets, eigen, kramtrx, params, error);
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
    // 使用tempk存储本征值
    VectorXcd tempk = eigen.k.segment(0, M);
    eigen.k.resize(M);
    eigen.k = tempk;
    // 使用tempphi存储本征函数
    MatrixXcd tempphi = eigen.phi.block(0, 0, eigen.M, eigen.phi.cols());
    eigen.phi.resize(M, eigen.phi.cols());
    eigen.phi = tempphi;
    // 使用tempphiR存储本征函数的实部
    MatrixXcd tempphiR = eigen.phiR.block(0, 0, eigen.M, eigen.phiR.cols());
    eigen.phiR.resize(M, eigen.phiR.cols());
    eigen.phiR = tempphiR;
    // 使用tempphiS存储本征函数的虚部
    MatrixXcd tempphiS = eigen.phiS.block(0, 0, eigen.M, eigen.phiS.cols());
    eigen.phiS.resize(M, eigen.phiS.cols());
    eigen.phiS = tempphiS;
    // 使用tempdphidz存储本征函数的导数
    MatrixXcd tempdphidz = eigen.dphidz.block(0, 0, eigen.M, eigen.dphidz.cols());
    eigen.dphidz.resize(M, eigen.dphidz.cols());
    eigen.dphidz = tempdphidz;
    // 使用tempdphidzR存储本征函数的导数的实部
    MatrixXcd tempdphidzR = eigen.dphidzR.block(0, 0, eigen.M, eigen.dphidzR.cols());
    eigen.dphidzR.resize(M, eigen.dphidzR.cols());
    eigen.dphidzR = tempdphidzR;
    // 使用tempdphidzS存储本征函数的导数的虚部
    MatrixXcd tempdphidzS = eigen.dphidzS.block(0, 0, eigen.M, eigen.dphidzS.cols());
    eigen.dphidzS.resize(M, eigen.dphidzS.cols());
    eigen.dphidzS = tempdphidzS;

    // // // 打印eigen.k
    // cout << "eigen.Extrap: \n"
    //      << eigen.Extrap.segment(0, M) << endl;

    // cout << "eigen.k: \n"
    //      << eigen.k << endl;

    for (int i = 0; i < M; i++)
    {
        eigen.k(i) = sqrt(eigen.Extrap(i) + eigen.k(i));
    }
}