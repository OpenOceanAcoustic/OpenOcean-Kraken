#ifndef SOLVE_H
#define SOLVE_H

#include "kkc_params.h"
#include "BCImpedanceMod.h"
#include "RootFinderBrent.h"
#include "RootFinderSecantMod.h"
#include "MergeVectors.h"
#include "InverseIteration.h"
#include "Scatter.h"

// 函数声明
void ERROUT();
void SolveEp(int &iset, const int &NSets, EigenParams &eigen, EigenFunction &eigenfun, KrakenMatrix &kramtrx, parameters &params, double &Error);
void Solve1(int &iset, const int& NSets, EigenParams &eigen, KrakenMatrix &kramtrx, parameters &params);
void Solve2(int &iset, EigenParams &eigen, KrakenMatrix &kramtrx, parameters &params);
void Solve3();
void FUNCT(int &iset, int &mode, double& x, double &Delta, int &iPower, KrakenMatrix &kramtrx,
           parameters& params, VectorXd &EVMat, bool coutmodes, int &modeCount);
void AcousticLayers(double x, double &f, double &g, int &iPower, KrakenMatrix &kramtrx, parameters& params, bool &coutmodes, int &modeCount);
void Bisection(int &iset, int &mode, double xMin, double xMax, VectorXd &xL, VectorXd &xR,
               KrakenMatrix &kramtrx, parameters &params, EigenParams &eigen);
void VectorSolve(KrakenMatrix &kramtrx, parameters &params, EigenParams &eigen, EigenFunction &eigenfun,
            int &NTotal, int &NTotal1);
void Normalize(int& mode, VectorXd &Phi, int &ITP, int &NTotal1, double &x, KrakenMatrix &kramtrx, parameters &params, EigenParams& eigen);
void ScatterLoss(int& mode, complex<double> &Perturbation_k, VectorXd &Phi, double &x, KrakenMatrix &kramtrx, parameters &params, EigenParams &eigen);

#endif // SOLVE_H