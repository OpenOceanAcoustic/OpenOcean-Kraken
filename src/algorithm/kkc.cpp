#include "kkc.h"

void run()
{
    parameters params;
    set_pekeris(params);

    // 初始化cp cs
    UpdateSSPLoss(params.freqinfo->freq, params.freqinfo->freq, params.NMedia, 
        params.SSPType, params.AttenUnit, params.SSP);
    UpdateHSLoss(params.freqinfo->freq, params.freqinfo->freq, params.NMedia, 
        params.AttenUnit, params.HSTop, params.HSBot);
    
    vector<int> NV{1, 2, 4, 8, 16};
    // TODO 计算本征值和本征函数
    KrakenMatrix kramtrx; 
    for (size_t iset = 0; iset < NSet; iset++)
    {
        int ntimes = NV[iset];
        Initialize(params, kramtrx, ntimes);
    }
}