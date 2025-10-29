#ifndef SOLVE_H
#define SOLVE_H

#include "kkc_params.h"
#include "BCImpedanceMod.h"

// 函数声明
void Solve1();
void FUNCT(int& mode, double x, double& Delta, int& iPower, KrakenMatrix& kramtrx, 
    parameters params, VectorXd& EVMat, bool& coutmodes, int& modeCount);
void AcousticLayers(double x, double& f, double& g, int& iPower);
void Bisection(double xMin, double xMax, std::vector<double>& xL, std::vector<double>& xR);

#endif // SOLVE_H