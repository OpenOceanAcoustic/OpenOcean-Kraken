#ifndef SOLVE_H
#define SOLVE_H

#include "kkc_params.h"

// 函数声明
void Solve1();
void FUNCT(double x, double& Delta, int& iPower);
void AcousticLayers(double x, double& f, double& g, int& iPower);
void Bisection(double xMin, double xMax, std::vector<double>& xL, std::vector<double>& xR);

// 外部变量声明（实际使用时可能需要调整）
extern int M, modeCount, ISet, FirstAcoustic, LastAcoustic, ifreq, iProf;
extern double omega2, cHigh, cLow;
extern std::vector<int> N, Loc;
extern std::vector<double> B1, B2, B3, B4, h, rho;
extern std::vector<std::vector<double>> EVMat;
extern std::vector<double> Extrap, k, VG;

// 结构体声明（可能需要根据项目实际情况调整）
struct SSPType {
    int NMedia;
};
extern SSPType* SSP;

#endif // SOLVE_H