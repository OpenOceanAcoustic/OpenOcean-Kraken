#include "kkc_params.h"
#include "src/set_pekeris.h"
#include "src/algorithm/init.h"
#include "src/algorithm/sspMod.h"

int main()
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
    Initialize(params, kramtrx);
    return 0;
}