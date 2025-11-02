#include "run.h"

void run()
{
    parameters params;
    set_pekeris(params);

    // 初始化cp cs
    UpdateSSPLoss(params.freqinfo->freq, params.freqinfo->freq, params.NMedia, 
        params.SSPType, params.AttenUnit, params.SSP);
    UpdateHSLoss(params.freqinfo->freq, params.freqinfo->freq, params.NMedia, 
        params.AttenUnit, params.HSTop, params.HSBot);
    
    // TODO 计算本征值和本征函数
    KrakenMatrix kramtrx; 
    EigenParams eigen;
    int NVsize = sizeof(params.mesh.NV)/sizeof(params.mesh.NV[0]);
    params.mesh.hV.resize(NVsize);
    for (int iset = 0; iset < NSet; iset++)
    {
        int ntimes = params.mesh.NV[iset];
        Initialize(iset, params, kramtrx, ntimes);
        Solve1(iset, NVsize, eigen, kramtrx, params);
        std::cout << "iset: " << iset << " \n" << eigen.EVMat.segment(iset*eigen.M, eigen.M).transpose() << std::endl;
    }
}