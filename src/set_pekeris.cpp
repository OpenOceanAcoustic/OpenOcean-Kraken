#include "set_pekeris.h"

void set_pekeris(parameters& params)
{
    params.freqinfo = new FreqInfo();
    params.Pos = new Position();

    // 介质层数
    params.NMedia = 1;
    params.SSP = new SSPStructure[params.NMedia];
    params.Bdry = new BdryType();

    // 频率
    params.freqinfo->freq = 100;

    // 海面、海底参数
    params.HSTop.BC = BC_Mode::MODE_V_Vacuum;
    params.HSBot.BC = BC_Mode::MODE_A_Half_space;
    params.HSBot.alphaI = 0.2;
    params.HSBot.alphaR = 1600.0;
    params.HSBot.betaI = 0.0;
    params.HSBot.betaR = 0.0;
    params.HSBot.Depth = 200.0;
    params.HSBot.rho = 1.5;
    params.HSBot.sigma = 0;

    // 相速度范围
    params.Chigh = 20000;
    params.Clow = 0;

    // 绝热模式
    params.modeType = ModeType::Adiabatic;

    // 声源接收设置
    params.Pos->NSz = 1;
    params.Pos->Sz.resize(params.Pos->NSz);
    params.Pos->Sz(0) = 25;
    params.Pos->NRr = 300;
    params.Pos->Rr.resize(params.Pos->NRr);
    for(size_t i = 0; i < params.Pos->NRr; ++i){
        params.Pos->Rr(i) = 100 * (i + 1);
    }
    params.Pos->NRz = 200;
    params.Pos->NRz_per_range = params.Pos->NRz;
    params.Pos->Rz.resize(params.Pos->NRz);
    params.Pos->Ro.resize(params.Pos->NRz);
    for(size_t i = 0; i < params.Pos->NRz; ++i){
        params.Pos->Rz(i) = i+1;
        params.Pos->Ro(i) = 0;
    }
    params.Pos->GridType = Grid_Mode::MODE_R_Rectangular;

    // 最大距离
    params.Rmax = 200000;

    // 计算模式，本征值和声场
    params.runMode = Run_Mode::MODE_B_Both;

    // 相干和非相干
    params.coherenceType = CoherenceType::Coherent;

    // 点声源
    params.SourceType = Source_Mode::MODE_R_Point;
    params.AttenUnit = Atten_Mode::MODE_W_db_per_lambda;

    // 声速剖面类型
    params.SSPType = SSP_Mode::MODE_C_cLinear;

    // 声速剖面
    params.SSP->NPts = 2;
    params.SSP->alphaR = Vector2d(1500.0, 1500.0);
    params.SSP->alphaI = Vector2d(0.0, 0.0);
    params.SSP->betaR = Vector2d(0.0, 0.0);
    params.SSP->betaI = Vector2d(0.0, 0.0);
    params.SSP->rho = Vector2d(1.0, 1.0);
    params.SSP->z = Vector2d(0, 200.0);
    params.SSP->cp.resize(params.SSP->NPts);
    params.SSP->cs.resize(params.SSP->NPts);
    params.SSP->beta = 0.0;
    params.SSP->ft = 0.0;
    params.SSP->N = 200;
    params.SSP->depth = 200;
    params.SSP->h = params.SSP->depth / (params.SSP->N - 1);
    params.SSP->sigma = 0;


    // params.SSP[1].NPts = 2;
    // params.SSP[1].alphaR = Vector2d(1600.0, 1600.0);
    // params.SSP[1].alphaI = Vector2d(0.1, 0.1);
    // params.SSP[1].betaR = Vector2d(0.0, 0.0);
    // params.SSP[1].betaI = Vector2d(0.0, 0.0);
    // params.SSP[1].rho = Vector2d(1.6, 1.6);
    // params.SSP[1].z = Vector2d(100, 120.0);
    // params.SSP[1].cp.resize(params.SSP[1].NPts);
    // params.SSP[1].cs.resize(params.SSP[1].NPts);
    // params.SSP[1].beta = 0.0;
    // params.SSP[1].ft = 1000;
    // params.SSP[1].N = 21;
    // params.SSP[1].depth = 20;
    // params.SSP[1].h = params.SSP[1].depth / (params.SSP[1].N - 1);

    params.Title = "Pekeris";
}
