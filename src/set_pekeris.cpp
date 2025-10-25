#include "set_pekeris.h"

void set_pekeris(parameters& params)
{
    params.freqinfo = new FreqInfo();
    params.Pos = new Position();
    params.SSP = new SSPStructure();
    params.Bdry = new BdryType();


    // 海面、海底参数
    params.HSTop.BC = BC_Mode::MODE_V_Vacuum;
    params.HSBot.BC = BC_Mode::MODE_A_Half_space;
    params.HSBot.alphaI = 0.2;
    params.HSBot.alphaR = 1600.0;
    params.HSBot.betaI = 0.0;
    params.HSBot.betaR = 0.0;
    params.HSBot.Depth = 100.0;
    params.HSBot.rho = 1.6;

    // 相速度范围
    params.Chigh = 2000;
    params.Clow = 1200;

    // 绝热模式
    params.modeType = ModeType::Adiabatic;

    // 介质层数
    params.NMedia = 1;

    // 声源接收设置
    params.Pos->NSz = 1;
    params.Pos->Sz.resize(params.Pos->NSz);
    params.Pos->Sz(0) = 25;
    params.Pos->NRr = 200;
    params.Pos->Rr.resize(params.Pos->NRr);
    for(size_t i = 0; i < params.Pos->NRr; ++i){
        params.Pos->Rr(i) = 100 * (i + 1);
    }
    params.Pos->NRz = 51;
    params.Pos->Rz.resize(params.Pos->NRr);
    for(size_t i = 0; i < params.Pos->NRr; ++i){
        params.Pos->Rz(i) = i*2;
    }

    // 最大距离
    params.Rmax = 0.0;

    // 计算模式，本征值和声场
    params.runMode = Run_Mode::MODE_B_Both;

    // 点声源
    params.SourceType = Source_Mode::MODE_R_Point;

    // 声速剖面
    params.SSP->Nz = 2;
    params.SSP->alphaR = Vector2d(1500.0, 1500.0);
    params.SSP->alphaI = Vector2d(0.0, 0.0);
    params.SSP->betaR = Vector2d(0.0, 0.0);
    params.SSP->betaI = Vector2d(0.0, 0.0);
    params.SSP->rho = Vector2d(1.0, 1.0);
    params.SSP->z = Vector2d(0, 100.0);
    params.SSP->AttenUnit = Atten_Mode::MODE_W_db_per_lambda;

    params.Title = "Pekeris";
}
