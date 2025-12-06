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
    params.freqinfo->freq = 150;
    params.freqinfo->Nfreq = 1;
    params.freqinfo->freqvec = VectorXd(1);
    params.freqinfo->freqvec(0) = params.freqinfo->freq;

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

void set_Munk(parameters& params)
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
    params.HSBot.alphaI = 0.5;
    params.HSBot.alphaR = 1506.5;
    params.HSBot.betaI = 0.0;
    params.HSBot.betaR = 0.0;
    params.HSBot.Depth = 3000.0;
    params.HSBot.rho = 1.5;
    params.HSBot.sigma = 0;

    // 相速度范围
    params.Chigh = 1600;
    params.Clow = 1400;

    // 绝热模式
    params.modeType = ModeType::Adiabatic;

    // 声源接收设置
    params.Pos->NSz = 1;
    params.Pos->Sz.resize(params.Pos->NSz);
    params.Pos->Sz(0) = 18.0;
    params.Pos->NRr = 1000;
    params.Pos->Rr.resize(params.Pos->NRr);
    for(int i = 0; i < params.Pos->NRr; ++i){
        params.Pos->Rr(i) = 100 * (i + 1);
    }
    params.Pos->NRz = 500;
    params.Pos->NRz_per_range = params.Pos->NRz;
    params.Pos->Rz.resize(params.Pos->NRz);
    params.Pos->Ro.resize(params.Pos->NRz);
    for(int i = 0; i < params.Pos->NRz; ++i){
        params.Pos->Rz(i) = 6.0 * (i+1);
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
    params.SSP->NPts = 23;
    params.SSP->alphaR.resize(23);
    params.SSP->alphaI = VectorXd::Ones(23)*0.5;
    params.SSP->betaR = VectorXd::Zero(23);
    params.SSP->betaI = VectorXd::Zero(23);
    params.SSP->rho = VectorXd::Ones(23);  // 水介质密度设为1.0
    params.SSP->z.resize(23);
    params.SSP->cp.resize(params.SSP->NPts);
    params.SSP->cs.resize(params.SSP->NPts);
    params.SSP->beta = 0.0;
    params.SSP->ft = 0.0;
    params.SSP->N = 200;  // 插值点数保持不变
    params.SSP->depth = 3000.0;  // 最大深度
    params.SSP->h = params.SSP->depth / (params.SSP->N - 1);
    params.SSP->sigma = 0;

    // 填充.env中的声速剖面数据
    params.SSP->z(0) = 0.0;     params.SSP->alphaR(0) = 1476.7;
    params.SSP->z(1) = 38.0;    params.SSP->alphaR(1) = 1476.7;
    params.SSP->z(2) = 50.0;    params.SSP->alphaR(2) = 1472.6;
    params.SSP->z(3) = 70.0;    params.SSP->alphaR(3) = 1468.8;
    params.SSP->z(4) = 100.0;   params.SSP->alphaR(4) = 1467.2;
    params.SSP->z(5) = 140.0;   params.SSP->alphaR(5) = 1471.6;
    params.SSP->z(6) = 160.0;   params.SSP->alphaR(6) = 1473.6;
    params.SSP->z(7) = 170.0;   params.SSP->alphaR(7) = 1473.6;
    params.SSP->z(8) = 200.0;   params.SSP->alphaR(8) = 1472.7;
    params.SSP->z(9) = 215.0;   params.SSP->alphaR(9) = 1472.2;
    params.SSP->z(10) = 250.0;  params.SSP->alphaR(10) = 1471.6;
    params.SSP->z(11) = 300.0;  params.SSP->alphaR(11) = 1471.6;
    params.SSP->z(12) = 370.0;  params.SSP->alphaR(12) = 1472.0;
    params.SSP->z(13) = 450.0;  params.SSP->alphaR(13) = 1472.7;
    params.SSP->z(14) = 500.0;  params.SSP->alphaR(14) = 1473.1;
    params.SSP->z(15) = 700.0;  params.SSP->alphaR(15) = 1474.9;
    params.SSP->z(16) = 900.0;  params.SSP->alphaR(16) = 1477.0;
    params.SSP->z(17) = 1000.0; params.SSP->alphaR(17) = 1478.1;
    params.SSP->z(18) = 1250.0; params.SSP->alphaR(18) = 1480.7;
    params.SSP->z(19) = 1500.0; params.SSP->alphaR(19) = 1483.8;
    params.SSP->z(20) = 2000.0; params.SSP->alphaR(20) = 1490.5;
    params.SSP->z(21) = 2500.0; params.SSP->alphaR(21) = 1498.3;
    params.SSP->z(22) = 3000.0; params.SSP->alphaR(22) = 1506.5;

    params.Title = "Dickins seamount";  
}