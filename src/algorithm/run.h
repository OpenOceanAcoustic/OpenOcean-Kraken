#ifndef RUN_H
#define RUN_H

#include "OpenOceanKrakenParams.h"
#include "ThreadPool.h"

#include "field.h"
#include "BCImpedanceMod.h"
#include "RootFinderMod.h"

#include "MergeVectors.h"
#include "InverseIteration.h"
#include "Scatter.h"
#include "sspMod.h"
#include <iomanip>

namespace OpenOceanKraken
{
    void EigenVWorker(ThreadPool &threadPool, const int &NumThreads, const size_t &iprof, const OOK_parameters &params, TridMtx &trid, OOK_output &output);
    void FieldWorker(const size_t &iprof, const OOK_parameters &params, OOK_output &output);

    // 函数声明

    void SolveEp(ThreadPool &threadPool, const int &NumThreads, const int &iset, const size_t &iprof, const int &NSets, EigenParams &eigen, const EigenParams *previousEigen, TridMtx &trid, const OOK_parameters &params, double &Error);
    void Solve1(ThreadPool &threadPool, const int &NumThreads, const int &iset, const size_t &iprof, const int &NSets, EigenParams &eigen, TridMtx &trid, const OOK_parameters &params);
    void Solve2(ThreadPool &threadPool, const int &NumThreads, const int &iset, const size_t &iprof, EigenParams &eigen, TridMtx &trid, const OOK_parameters &params);
    bool Solve3(ThreadPool &threadPool, const int &NumThreads, const int &iset, const size_t &iprof, EigenParams &eigen, const EigenParams &previousEigen, TridMtx &trid, const OOK_parameters &params);
    void TridPreprocess(int &iset, size_t iprof, const OOK_parameters &params, TridMtx &trid, int ntimes);

    void VectorSolve(size_t iprof, TridMtx &trid, const OOK_parameters &params, EigenParams &eigen,
                     int &NTotal, int &NTotal1);
    void Normalize(const size_t &iprof, const int &mode, const int &firstM, Eigen::VectorXd &Phi, int &ITP, int &NTotal1, double &x, TridMtx &trid, const OOK_parameters &params, EigenParams &eigen);
    void ScatterLoss(const size_t &iprof, const int &mode, std::complex<double> &Perturbation_k, Eigen::VectorXd &Phi, double &x, TridMtx &trid, const OOK_parameters &params, EigenParams &eigen);

}
#endif
