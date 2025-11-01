#ifndef SOLVE_H
#define SOLVE_H

#include "kkc_params.h"
#include "BCImpedanceMod.h"
#include "RootFinderBrent.h"

// 函数声明
void ERROUT(const std::string &routine, const std::string &message);

void Solve1(int &iset, const int& NSets, EigenParams &eigen, KrakenMatrix &kramtrx, parameters &params);
void FUNCT(int &iset, int &mode, double& x, double &Delta, int &iPower, KrakenMatrix &kramtrx,
           parameters& params, VectorXd &EVMat, bool coutmodes, int &modeCount);
void AcousticLayers(double x, double &f, double &g, int &iPower, KrakenMatrix &kramtrx, bool &coutmodes, int &modeCount);
void Bisection(int &iset, int &mode, double xMin, double xMax, VectorXd &xL, VectorXd &xR,
               KrakenMatrix &kramtrx, parameters &params, EigenParams &eigen);

#endif // SOLVE_H