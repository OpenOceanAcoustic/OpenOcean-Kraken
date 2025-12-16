#ifndef SOLVE_H
#define SOLVE_H

#include "kkc_params.h"
#include "BCImpedanceMod.h"
#include "RootFinderBrent.h"
#include "RootFinderSecantMod.h"
#include "MergeVectors.h"
#include "InverseIteration.h"
#include "Scatter.h"
#include <iomanip>

// 函数声明
void ERROUT();
void SolveEp(int &iset, size_t iprof, const int &NSets, EigenParams &eigen, EigenFunction &eigenfun, KrakenMatrix &kramtrx, parameters &params, double &Error);
void Solve1(int &iset, size_t iprof, const int& NSets, EigenParams &eigen, KrakenMatrix &kramtrx, parameters &params);
void Solve2(int &iset, size_t iprof, EigenParams &eigen, KrakenMatrix &kramtrx, parameters &params);
void Solve3(int &iset, size_t iprof, EigenParams &eigen, KrakenMatrix &kramtrx, parameters &params);
void FUNCT(int &iset, size_t iprof, int &mode, double& x, double &Delta, int &iPower, KrakenMatrix &kramtrx,
           parameters& params, VectorXd &EVMat, const int& firstM, bool coutmodes, int &modeCount);
void AcousticLayers(double x, double &f, double &g, int &iPower, KrakenMatrix &kramtrx, parameters& params, bool &coutmodes, int &modeCount);
void Bisection(int &iset, size_t iprof, int &mode, double xMin, double xMax, VectorXd &xL, VectorXd &xR,
               KrakenMatrix &kramtrx, parameters &params, EigenParams &eigen);
void VectorSolve(size_t iprof, KrakenMatrix &kramtrx, parameters &params, EigenParams &eigen, EigenFunction &eigenfun,
            int &NTotal, int &NTotal1);
void Normalize(size_t iprof, int& mode, int& firstM, VectorXd &Phi, int &ITP, int &NTotal1, double &x, KrakenMatrix &kramtrx, parameters &params, EigenParams& eigen);
void ScatterLoss(size_t iprof, int& mode, complex<double> &Perturbation_k, VectorXd &Phi, double &x, KrakenMatrix &kramtrx, parameters &params, EigenParams &eigen);

#endif // SOLVE_H