#ifndef SOLVE_H
#define SOLVE_H

#include "kkc_params.h"
#include "BCImpedanceMod.h"
#include "RootFinderBrent.h"
#include "RootFinderSecantMod.h"
#include "MergeVectors.h"
#include "InverseIteration.h"
#include "Scatter.h"
#include "sspMod.h"
#include <iomanip>

// 函数声明
void ERROUT();
void SolveEp(const int &iset, const size_t& iprof, const int &NSets, EigenParams &eigen, TridMtx &trid, parameters &params, double &Error);
void Solve1(const int &iset, const size_t& iprof, const int& NSets, EigenParams &eigen, TridMtx &trid, parameters &params);
void Solve2(const int &iset, const size_t& iprof, EigenParams &eigen, TridMtx &trid, parameters &params);
void Solve3(const int &iset, const size_t& iprof, EigenParams &eigen, TridMtx &trid, parameters &params);
void TridPreprocess(int &iset, size_t iprof, parameters &params, TridMtx &trid, int ntimes);
void FUNCT(const int &iset, const size_t& iprof, int &mode, double& x, double &Delta, int &iPower, TridMtx &trid,
           parameters& params, VectorXd &EVMat, const int& firstM, const bool &isCountMode, int &modeCount);
void AcousticLayers(const size_t& iprof, double x, double &f, double &g, int &iPower, TridMtx &trid, parameters& params, const bool &isCountMode, int &modeCount);
void Bisection(const int &iset, const size_t& iprof, int &mode, double xMin, double xMax, VectorXd &xL, VectorXd &xR,
               TridMtx &trid, parameters &params, EigenParams &eigen);
void VectorSolve(size_t iprof, TridMtx &trid, parameters &params, EigenParams &eigen,
                 int &NTotal, int &NTotal1);
void Normalize(const size_t& iprof, const int& mode, const int& firstM, VectorXd &Phi, int &ITP, int &NTotal1, double &x, TridMtx &trid, parameters &params, EigenParams& eigen);
void ScatterLoss(const size_t& iprof, const int& mode, complex<double> &Perturbation_k, VectorXd &Phi, double &x, TridMtx &trid, parameters &params, EigenParams &eigen);

#endif // SOLVE_H