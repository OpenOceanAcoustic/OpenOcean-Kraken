#include "run.h"
#include <cmath>
#include <limits>
#include <stdexcept>

// 计算本征值和本征函数
namespace OpenOceanKraken
{
    namespace
    {
        int MinLocMaskedPositive(const Eigen::VectorXd &values, const int start, const int count, const double threshold)
        {
            int minLoc = 0;
            double minValue = std::numeric_limits<double>::infinity();
            for (int mode = 0; mode < count; ++mode)
            {
                const double value = values(start + mode);
                if (std::isfinite(value) && value > threshold && value < minValue)
                {
                    minValue = value;
                    minLoc = mode + 1;
                }
            }
            return minLoc;
        }

        std::complex<double> HalfSpaceRadiationRoot(const std::complex<double> &z)
        {
            std::complex<double> root = std::sqrt(z);
            if (std::imag(root) < 0.0)
            {
                root = -root;
            }
            return root;
        }
    }

    void EigenVWorker(ThreadPool &threadPool, const int &NumThreads, const size_t &iprof, const OOK_parameters &params, TridMtx &trid, OOK_output &output)
    {
        double freq = params.freqinfo.freq;
        auto &ssp = params.SSP.at(iprof);
        auto &HSTop = params.SSP.at(iprof).HSTop;
        auto &HSBot = params.SSP.at(iprof).HSBot;
        auto &eigen = output.eigen[iprof];
        const EigenParams *previousEigen = nullptr;
        if (iprof > 0 && params.modeType == ModeType::Couple)
        {
            previousEigen = &output.eigen[iprof - 1];
        }
        double error;

        for (int iset = 0; iset < params.mesh.NSets; iset++)
        {
            int ntimes = params.mesh.NV[iset];
            TridPreprocess(iset, iprof, params, trid, ntimes);
            SolveEp(threadPool, NumThreads, iset, iprof, params.mesh.NSets, eigen, previousEigen, trid, params, error);
            if (error * params.Rmax < 1.0)
            {
                break;
            }
            else
            {
                if (iset == params.mesh.NSets - 1)
                    std::cout << "Warning in KRAKEN : Too many meshes needed: check convergence" << std::endl;
            }
        }

        const double finalThreshold = SQ(2 * pi * params.freqinfo.freq / trid.cHigh);
        eigen.M = MinLocMaskedPositive(eigen.Extrap, 0, eigen.M, finalThreshold);
        if (eigen.M <= 0)
        {
            throw std::runtime_error("KRAKEN eigen solve produced zero propagating modes.");
        }

        for (int i = 0; i < eigen.M; i++)
        {
            eigen.k(i) = sqrt(eigen.Extrap(i) + eigen.k(i));
        }
    }

    // 计算声场
    void FieldWorker(const size_t &iprof, const OOK_parameters &params, OOK_output &output)
    {
        auto eigen = output.eigen[iprof];
        eigen.M = std::min(eigen.M, params.MLimit);
        for (int isz = 0; isz < params.Pos.NSz; isz++)
        {
            Evaluate(eigen, params, isz, iprof, output.u_AllSources, output.v_AllSources, output.h_AllSources);
        }
    }

    void SolveEp(ThreadPool &threadPool, const int &NumThreads, const int &iset, const size_t &iprof, const int &NSets, EigenParams &eigen, const EigenParams *previousEigen, TridMtx &trid, const OOK_parameters &params, double &Error)
    {
        double omega2 = SQ(2 * pi * params.freqinfo.freq);
        const auto &ssp_prof = params.SSP.at(iprof); // 只查一次
        if (iprof > 0 && iset < 2 && params.modeType == ModeType::Couple)
        {
            if (previousEigen == nullptr)
            {
                throw std::logic_error("KRAKEN Solve3 requires the previous profile eigenvalues.");
            }
            if (!Solve3(threadPool, NumThreads, iset, iprof, eigen, *previousEigen, trid, params))
            {
                if (ssp_prof.NMedia <= ssp_prof.LastAcoustic - ssp_prof.FirstAcoustic + 1)
                {
                    std::cout << "Warning in KRAKEN - Solve3: previous profile initial guesses are unavailable; falling back to Solve1." << std::endl;
                    Solve1(threadPool, NumThreads, iset, iprof, NSets, eigen, trid, params);
                }
                else
                {
                    throw std::runtime_error("KRAKEN Solve3 has insufficient valid initial guesses for an elastic profile.");
                }
            }
        }
        else if ((iset < 2) && (ssp_prof.NMedia <= ssp_prof.LastAcoustic - ssp_prof.FirstAcoustic + 1))
        {
            Solve1(threadPool, NumThreads, iset, iprof, NSets, eigen, trid, params);
        }
        else
        {
            Solve2(threadPool, NumThreads, iset, iprof, eigen, trid, params);
        }

        int start_idx = iset * eigen.firstM;
        eigen.Extrap.segment(start_idx, eigen.M) = eigen.EVMat.segment(start_idx, eigen.M);

        // std::cout << "iset: " << iset << " \n"
        //           << eigen.EVMat.segment(start_idx, eigen.M) << std::endl;

        // 查找满足条件的最小位置
        double threshold = omega2 / SQ(trid.cHigh);

        eigen.M = MinLocMaskedPositive(eigen.Extrap, 0, eigen.M, threshold);

        // 计算NTotal：N(FirstAcoustic : LastAcoustic)的和
        int NTotal = 0;
        for (int i = ssp_prof.FirstAcoustic; i <= ssp_prof.LastAcoustic; ++i)
        {
            NTotal += trid.N(i);
        }
        int NTotal1 = NTotal + 1;

        // 如果是第一个网格，计算特征向量
        if (iset == 0)
        {
            VectorSolve(iprof, trid, params, eigen, NTotal, NTotal1);
        }

        // // 打印特征值
        // std::cout << "iset: " << iset << " \n" << std::endl;
        // for (int mode = 0; mode < eigen.M; ++mode)
        // {
        //     // 8位精度输出
        //     std::cout << std::setprecision(8) << eigen.EVMat(iset * eigen.firstM + mode) << std::endl;
        // }

        // 初始化误差和KEY
        Error = 1.0e10;
        int KEY = 2 * eigen.M / 3;

        // 1. x2 只依赖 iset，与内层循环无关 → 提到最外层
        const double x2_sq = SQ(params.mesh.NV[iset]);
        if (iset > 0)
        {
            // T1 = Extrap(1, KEY)，对应索引KEY-1
            double T1 = eigen.Extrap(KEY);

            // 理查森外推改善精度
            for (int j = iset - 1; j >= 0; --j)
            {
                const double x1_sq = SQ(params.mesh.NV[j]); // ← 每个 j 算一次
                const double ratio = x2_sq / x1_sq;         // 避免重复除法
                const double denom = ratio - 1.0;
                for (int mode = 0; mode < eigen.M; ++mode)
                {
                    // 计算当前(j, mode)和(j+1, mode)在向量中的索引
                    int idx_j = j * eigen.firstM + mode;
                    int idx_j1 = (j + 1) * eigen.firstM + mode;

                    double F1 = eigen.Extrap(idx_j);
                    double F2 = eigen.Extrap(idx_j1);

                    // 理查森外推公式
                    // eigen.Extrap(idx_j) = F2 - (F1 - F2) / (x2 / x1 - 1.0);
                    eigen.Extrap(idx_j) = F2 - (F1 - F2) / denom;
                }
            }

            // 计算误差
            double T2 = eigen.Extrap(KEY);
            Error = abs(T2 - T1);
        }
    }

    // Solve1函数：使用Sturm序列分离本征值以及用Brent求根法求得本征值
    void Solve1(ThreadPool &threadPool, const int &NumThreads, const int &iset, const size_t &iprof, const int &NSets, EigenParams &eigen, TridMtx &trid, const OOK_parameters &params)
    {
        int iPower = 0, NTotal, mode = 0; // NzTab = 0,
        double x, x1, x2, xMin, xMax, Eps, Delta;
        bool isCountMode = true;
        std::string ErrorMessage;
        Eigen::VectorXd xL, xR;
        double omega2 = SQ(2 * pi * params.freqinfo.freq);

        int modeCount = 0;

        // 确定模态数量
        xMin = 1.00001 * omega2 / SQ(trid.cHigh);

        FUNCT(iset, iprof, mode, xMin, Delta, iPower, trid, params, eigen.EVMat, eigen.firstM, isCountMode, modeCount);
        int M = modeCount;
        eigen.M = M;

        // 分配xL和xR的内存
        xL.resize(M + 1);
        xR.resize(M + 1);

        xMax = omega2 / SQ(trid.cLow); // 最大波数的平方
        FUNCT(iset, iprof, mode, xMax, Delta, iPower, trid, params, eigen.EVMat, eigen.firstM, isCountMode, modeCount);

        M = M - modeCount;
        eigen.M = M;

        if (M == 0)
        {
            // 处理无模态情况
            // LRecordLength = 32; // MOD文件片段长度

            // 打开MODFile并写入头部
            // 在C++中需要使用适当的文件I/O

            // 抛出错误
        }

        // 计算NTotal
        NTotal = 0;
        for (int i = params.SSP[iprof].FirstAcoustic; i <= params.SSP[iprof].LastAcoustic; ++i)
        {
            NTotal += trid.N(i);
        }

        if (M > NTotal / 5)
        {
            // 输出警告信息
            std::cout << "Approximate number of modes = " << M << " Warning in KRAKEN - Solve1 : Mesh too coarse to sample the modes adequately" << std::endl;
        }

        Bisection(iset, iprof, mode, xMin, xMax, xL, xR, trid, params, eigen); // 初始化上下边界

        // 使用ZBRENT精化每个本征值
        isCountMode = false;
        // 缓存求解结果（根值和是否成功）
        std::vector<double> roots(M, 0.0);
        std::vector<unsigned char> rootOK(M, 0);

        // 并行处理每个模态
        for (int threadId = 0; threadId < NumThreads; ++threadId)
        {
            threadPool.enqueue_with_id("RunZBRENTX", [&NumThreads, &params, &eigen, &trid, &xL, &xR, &roots, &rootOK, M, iset, iprof, mode, isCountMode, threadId]()
                                       {
            // 当前线程处理的任务 ID 满足 modeIdx % NumThreads == threadId
            for (int modeIdx = threadId; modeIdx < M; modeIdx += NumThreads)
            {
                double local_x, local_x1, local_x2, local_Eps, local_Delta;
                int local_iPower = 0, local_modeCount = 0;

                local_x1 = xL(modeIdx);
                local_x2 = xR(modeIdx);
                local_x = local_x1;
                local_Eps = std::abs(local_x2) * std::pow(10.0, 2.0 - std::numeric_limits<double>::digits10);
                
                bool ok = ZBRENTX(local_x, local_x1, local_x2, local_Eps, iset, iprof, mode, local_Delta, local_iPower, trid, params, eigen.EVMat, eigen.firstM, isCountMode, local_modeCount);

                // 直接判断 ok，不再判断 local_ErrorMessage
                if (!ok)
                {
                    // 构造一个包含上下文的警告信息，因为具体的错误原因拿不到
                    std::string warningMsg = "ZBRENTX 求解失败 (模态索引: " + std::to_string(modeIdx) +
                                             ", 区间: [" + std::to_string(local_x1) +
                                             ", " + std::to_string(local_x2) + "])";

                    // 输出到日志或控制台
                    std::cerr << "[WARNING] " << warningMsg << std::endl;

                }
                else
                {
                    roots[modeIdx] = local_x;
                    rootOK[modeIdx] = 1;
                }
            } });
        }

        // 等待所有ZBRENTX任务完成
        threadPool.wait_id("RunZBRENTX");

        int failedModes = 0;
        for (int modeIdx = 0; modeIdx < M; ++modeIdx)
        {
            if (rootOK[modeIdx])
            {
                eigen.EVMat(iset * eigen.firstM + modeIdx) = roots[modeIdx];
            }
            else
            {
                ++failedModes;
            }
        }
        if (failedModes != 0)
        {
            std::cout << "Warning in KRAKEN - Solve1 : Brent failed for "
                      << failedModes << " isolated mode bracket(s)" << std::endl;
        }
    }

    void Solve2(ThreadPool &threadPool, const int &NumThreads, const int &iset, const size_t &iprof, EigenParams &eigen, TridMtx &trid, const OOK_parameters &params)
    {
        double omega2 = SQ(2 * pi * params.freqinfo.freq);
        double x1, x2, Tolerance, Delta;
        double x = omega2 / SQ(trid.cLow);
        int Iteration, MaxIteration = 2000, iPower = 0, modeCount = 0;

        bool isCountMode = false;

        Eigen::VectorXd P(10);
        if (eigen.firstM <= 0)
        {
            throw std::runtime_error("KRAKEN Solve2 requires a positive modal capacity.");
        }
        if (eigen.M <= 0 || eigen.M > eigen.firstM)
        {
            eigen.M = eigen.firstM;
        }

        for (int mode = 0; mode < eigen.M; mode++)
        {
            x = 1.00001 * x;
            if (iset >= 1)
            {
                for (int i = 0; i < iset; i++)
                {
                    P(i) = eigen.EVMat(i * eigen.firstM + mode);
                }
                if (iset >= 2)
                {
                    for (int ii = 0; ii < iset - 1; ii++)
                    {
                        for (int j = 0; j < iset - ii - 1; j++)
                        {
                            x1 = SQ(trid.hV(j));
                            x2 = SQ(trid.hV(j + ii + 1));
                            P(j) = ((SQ(trid.hV(iset)) - x2) * P(j) -
                                    (SQ(trid.hV(iset)) - x1) * P(j + 1)) /
                                   (x1 - x2);
                        }
                    }
                    x = P(0);
                }
            }
            Tolerance = abs(x) * trid.B1.size() * pow(10.0, (1.0 - std::numeric_limits<double>::digits10));
            Iteration = MaxIteration + 1;
            ZSecantX(x, Tolerance, Iteration, MaxIteration, iset, iprof, mode, Delta, iPower, trid, params, eigen.EVMat, eigen.firstM, isCountMode, modeCount);
            if (Iteration > MaxIteration || !std::isfinite(x))
            {
                std::cout << "Warning in KRAKEN - Solve2 : RootFinderSecant failed for mode "
                          << (mode + 1) << std::endl;
                x = std::numeric_limits<double>::min();
            }
            eigen.EVMat(iset * eigen.firstM + mode) = x;
            if (omega2 / SQ(trid.cHigh) > x)
            {
                eigen.M = mode;
                return;
            }
        }

    }

    bool Solve3(ThreadPool &threadPool, const int &NumThreads, const int &iset, const size_t &iprof, EigenParams &eigen, const EigenParams &previousEigen, TridMtx &trid, const OOK_parameters &params)
    {
        int IT, MaxIT, iPower = 0, mode = 0; // NzTab = 0,
        double x, xMin, Tolerance, Delta;

        bool isCountMode = true;

        MaxIT = 500;

        double omega2 = SQ(2 * pi * params.freqinfo.freq);

        int modeCount = 0;

        // 确定模态数量
        xMin = 1.00001 * omega2 / SQ(trid.cHigh);

        FUNCT(iset, iprof, mode, xMin, Delta, iPower, trid, params, eigen.EVMat, eigen.firstM, isCountMode, modeCount);
        int M = modeCount;
        eigen.M = M;
        isCountMode = false;

        const Eigen::Index previousSetStart = static_cast<Eigen::Index>(iset) * previousEigen.firstM;
        const Eigen::Index previousSetEnd = previousSetStart + M;
        if (M <= 0 || M > eigen.firstM || M > previousEigen.firstM || previousSetStart < 0 || previousSetEnd > previousEigen.EVMat.size())
        {
            return false;
        }

        for (int modeIdx = 0; modeIdx < M; ++modeIdx)
        {
            x = previousEigen.EVMat(previousSetStart + modeIdx);
            if (!std::isfinite(x) || x <= 0.0)
            {
                return false;
            }
            Tolerance = std::abs(x) * std::pow(10.0, 2.0 - std::numeric_limits<double>::digits10);
            IT = MaxIT + 1;
            ZSecantX(x, Tolerance, IT, MaxIT, iset, iprof, modeIdx, Delta, iPower, trid, params, eigen.EVMat, eigen.firstM, isCountMode, modeCount);

            if (IT > MaxIT || !std::isfinite(x))
            {
                std::cout << "Warning in KRAKEN - Solve3: RootFinderSecant failed for profile "
                          << (iprof + 1) << ", mesh set " << (iset + 1)
                          << ", mode " << (modeIdx + 1) << std::endl;
                x = std::numeric_limits<double>::min();
                // 输出警告信息
            }

            eigen.EVMat(iset * eigen.firstM + modeIdx) = x;

            if (omega2 / SQ(trid.cHigh) > x)
            {
                eigen.M = modeIdx; // 调整为当前索引
                return true;
            }
        }
        eigen.M = M;
        return true;
    }

    // 初始化有限差分方程
    void TridPreprocess(int &iset, size_t iprof, const OOK_parameters &params, TridMtx &trid, int ntimes)
    {
        bool ElasticFlag = false;
        int NPoints = 0;

        double Two_h;
        double cp2, cs2;
        double omega2 = SQ(2 * pi * params.freqinfo.freq);

        // 初始化变量
        double cLow = params.cLow;
        double cMin = 1e8;
        double cHigh = params.cHigh;
        trid.Loc[0] = 0; // C++使用0-based索引

        auto ssp_profile = params.SSP.at(iprof);
        auto HSTop = ssp_profile.HSTop;
        auto HSBot = ssp_profile.HSBot;

        for (int i = 0; i < ssp_profile.NMedia; ++i)
        {
            trid.N(i) = ssp_profile.NMesh(i) * ntimes;
            trid.h(i) = ssp_profile.depth(i) / trid.N(i);
            // 计算当前层的起始位置
            if (i != 0)
            {
                trid.Loc(i) = trid.Loc(i - 1) + trid.N(i - 1) + 1;
            }
            NPoints += trid.N(i);
            if (i == 0)
            {
                trid.hV(iset) = trid.h(i);
            }
            ssp_profile.NMesh(i) = trid.N(i);
            ssp_profile.interp_offset(i) = trid.Loc(i);
        }
        NPoints += ssp_profile.NMedia;

        // 处理每个介质层
        for (int im = 0; im < ssp_profile.NMedia; ++im)
        {
            int ii = trid.Loc(im); // C++使用0-based索引，不需要+1
            Two_h = 2.0 * trid.h(im);

            EvaluateSSP(trid, ssp_profile, ssp_profile.SSPType, im);

            // 打印rho_int
            // std::cout << "RHO INT" << SSP.rho_int.transpose() << std::endl;

            // 加载有限差分方程的对角线
            if (ssp_profile.Material.at(im) == Media_Mode::MODE_A_Acoustic)
            { // 声学介质情况

                // 计算当前层的最小声速
                double min_cp = 1e8;
                for (int j = 0; j < trid.N(im) + 1; ++j)
                {
                    min_cp = std::min(min_cp, std::real(trid.cp_int(ii + j)));
                }
                cMin = std::min(cMin, min_cp);

                // 计算B1和B1C
                double h2 = SQ(trid.h(im));
                for (int j = 0; j < trid.N(im) + 1; ++j)
                {
                    std::complex<double> val = omega2 / SQ(trid.cp_int(ii + j));
                    trid.B1(ii + j) = -2.0 + h2 * std::real(val);
                    trid.B1C(ii + j) = std::imag(val);
                    trid.rho(ii + j) = trid.rho_int(ii + j);
                }
            }
            else
            { // 弹性介质情况
                ElasticFlag = true;
                Two_h = 2.0 * trid.h(im);

                for (int j = 0; j < trid.N(im) + 1; ++j)
                {
                    const int idx = ii + j;
                    cMin = std::min(std::real(trid.cs_int(idx)), cMin);

                    cp2 = std::real(trid.cp_int(idx) * trid.cp_int(idx));
                    cs2 = std::real(trid.cs_int(idx) * trid.cs_int(idx));

                    trid.B1(idx) = Two_h / (trid.rho_int(idx) * cs2);
                    trid.B2(idx) = Two_h / (trid.rho_int(idx) * cp2);
                    trid.B3(idx) = 4.0 * Two_h * trid.rho_int(idx) * cs2 * (cp2 - cs2) / cp2;
                    trid.B4(idx) = Two_h * (cp2 - 2.0 * cs2) / cp2;
                    trid.rho(idx) = Two_h * std::real(omega2) * trid.rho_int(idx);
                }
            }
        }

        // // 打印kramtrx内部参数
        // std::cout << "B1 " << kramtrx.B1.size() << kramtrx.B1.transpose() << std::endl;

        // std::cout << "B1C " << kramtrx.B1C.size() << kramtrx.B1C.transpose() << std::endl;

        // std::cout << "rho " << kramtrx.rho.size() << kramtrx.rho.transpose() << std::endl;

        // 处理底部半空间属性
        if (HSBot.BC == BC_Mode::MODE_A_Half_space)
        {
            if (std::real(HSBot.cs) > 0.0)
            { // 弹性底部半空间
                ElasticFlag = true;
                cMin = std::min(cMin, std::real(HSBot.cs));
                cHigh = std::min(cHigh, std::real(HSBot.cs));
            }
            else
            { // 声学底部半空间
                cMin = std::min(cMin, std::real(HSBot.cp));
            }
        }

        // 处理顶部半空间属性
        if (HSTop.BC == BC_Mode::MODE_A_Half_space)
        {
            if (std::real(HSTop.cs) > 0.0)
            { // 弹性顶部半空间
                ElasticFlag = true;
                cMin = std::min(cMin, std::real(HSTop.cs));
                cHigh = std::min(cHigh, std::real(HSTop.cs));
            }
            else
            { // 声学顶部半空间
                cMin = std::min(cMin, std::real(HSTop.cp));
            }
        }

        // 如果存在弹性介质，则减小cMin以考虑Scholte波
        if (ElasticFlag)
        {
            cMin = 0.85 * cMin;
        }
        trid.cLow = std::max(cLow, cMin);
        trid.cHigh = cHigh;
    }

    void VectorSolve(size_t iprof, TridMtx &trid, const OOK_parameters &params, EigenParams &eigen,
                     int &NTotal, int &NTotal1)
    {
        // === 1. 缓存频繁访问的引用 ===
        const auto &ssp = params.SSP[iprof]; // 假设 iprof 合法；若需安全检查，用 .at(iprof)
        const int first_ac = ssp.FirstAcoustic;
        const int last_ac = ssp.LastAcoustic;
        const auto &z_ssp = ssp.z;

        // === 2. 构建深度网格 z 和系数 e ===
        Eigen::VectorXd z(NTotal1), e(NTotal1 + 1);
        int j = 0;
        double h_rho = 0.0;

        z(0) = z_ssp[first_ac];
        for (int im = first_ac; im <= last_ac; ++im)
        {
            const double h_im = trid.h(im);
            h_rho = h_im * trid.rho(trid.Loc(im));
            for (int ii = 1; ii <= trid.N(im); ++ii)
            {
                e(j + ii) = 1.0 / h_rho;
                z(j + ii) = z(j) + h_im * ii;
            }
            j += trid.N(im);
        }
        e(NTotal1) = 1.0 / h_rho;

        // === 3. 一次性构建插值表（关键优化：移出 mode 循环！）===
        const Position &modePos = params.hasModePos ? params.ModePos : params.Pos;
        Eigen::VectorXd zTab;
        int NzTab;
        Eigen::VectorXi Ix, Iy;
        MergeVectors(modePos.Sz, modePos.Rz, zTab, NzTab, Ix, Iy);

        Eigen::VectorXd WTS(params.Pos.NSz), WTR(params.Pos.NRz), WTZ(NzTab);
        Eigen::VectorXi ISzTab(params.Pos.NSz), IRzTab(params.Pos.NRz), IZzTab(NzTab);

        Weight_dble(z, NTotal1, zTab, NzTab, WTZ, IZzTab);
        Weight_dble(zTab, NzTab, params.Pos.Sz, params.Pos.NSz, WTS, ISzTab);
        Weight_dble(zTab, NzTab, params.Pos.Rz, params.Pos.NRz, WTR, IRzTab);

        // 使用 float 矩阵直接匹配 .mod 输出格式
        eigen.ModeZ = zTab;
        eigen.PhiMode.resize(eigen.firstM, NzTab);
        eigen.PhiMode.setZero();
        Eigen::MatrixXcd dPhiMode(eigen.firstM, NzTab);
        dPhiMode.setZero();

        // === 4. 临时工作向量（在 mode 循环内重用）===
        Eigen::VectorXd Psi(NTotal1), dPsidz(NTotal1), d(NTotal1);
        Eigen::VectorXd dPsi(NTotal1 - 1); // 差分缓存

        // === 5. 主循环：每个模式 ===
        for (int mode = 0; mode < eigen.M; ++mode)
        {
            double x = eigen.EVMat(mode);
            std::complex<double> fTop, gTop, fBot, gBot;
            bool isComplex = false;
            int iPower = 0, modeCount = 0;

            // --- 上边界条件 ---
            BCImpedance(iprof, x, true, fTop, gTop, iPower, isComplex, trid, params, modeCount);
            if (gTop == 0.0)
            {
                d(0) = 1.0;
                e(1) = 0.0;
            }
            else
            {
                int L = trid.Loc(first_ac);
                double xh2 = x * SQ(trid.h(first_ac));
                h_rho = trid.h(first_ac) * trid.rho(L);
                d(0) = (trid.B1(L) - xh2) / h_rho / 2.0 + std::real(fTop / gTop);
            }

            // --- 构建三对角矩阵主对角线 d ---
            int ITP = NTotal;
            j = 0;
            int L = trid.Loc(first_ac);
            for (int im = first_ac; im <= last_ac; ++im)
            {
                const double h_im = trid.h(im);
                const double xh2 = x * SQ(h_im);
                // 注意：rho 索引从 L+1 开始（因 z 从界面开始）
                h_rho = h_im * trid.rho(L);

                if (im > first_ac)
                {
                    L++;
                    d(j) = (d(j) + (trid.B1(L) - xh2) / h_rho) / 2.0;
                }

                for (int ii = 0; ii < trid.N(im); ++ii)
                {
                    j++;
                    L++;
                    d(j) = (trid.B1(L) - xh2) / h_rho;
                    if (trid.B1(L) - xh2 + 2.0 > 0.0)
                    {
                        ITP = std::min(j, ITP);
                    }
                }
            }

            // --- 下边界条件 ---
            BCImpedance(iprof, x, false, fBot, gBot, iPower, isComplex, trid, params, modeCount);
            if (gBot == 0.0)
            {
                d(NTotal1 - 1) = 1.0;
                e(NTotal1 - 1) = 0.0;
            }
            else
            {
                d(NTotal1 - 1) = d(NTotal1 - 1) / 2.0 - std::real(fBot / gBot);
            }

            // --- 求解特征向量 ---
            int IErr = 0;
            InverseIterationD(NTotal1, d, e, IErr, Psi);

            if (IErr != 0)
            {
#ifndef NDEBUG
                std::cerr << "Warning in KRAKEN - InverseIteration: mode=" << mode
                          << " failed to converge.\n";
#endif
                Psi.setZero(); // 清零失败模式
            }
            else
            {
                Normalize(iprof, mode, eigen.firstM, Psi, ITP, NTotal1, x, trid, params, eigen);
            }

            // --- 计算 dPsi/dz（使用中心差分）---
            // 预计算差分 Psi[i+1] - Psi[i]
            for (int idx = 0; idx < NTotal1 - 1; ++idx)
            {
                dPsi(idx) = Psi(idx + 1) - Psi(idx);
            }

            // 边界点
            dPsidz(0) = dPsi(0) / trid.h(first_ac);
            dPsidz(NTotal) = dPsi(NTotal - 1) / trid.h(last_ac);

            // 内部点（注意：需根据深度确定当前 h）
            int idx_z = 1; // 当前全局 z 索引
            for (int im = first_ac; im <= last_ac; ++im)
            {
                const double h_im = trid.h(im);
                const int N_seg = (im == last_ac) ? trid.N(im) + 1 : trid.N(im);
                for (int ii = 0; ii < N_seg; ++ii)
                {
                    if (idx_z > 0 && idx_z < NTotal)
                    {
                        // 中心差分：需要知道左右 h，简化为使用当前层 h
                        dPsidz(idx_z) = dPsi(idx_z - 1) / (2.0 * h_im) + dPsi(idx_z) / (2.0 * h_im);
                    }
                    idx_z++;
                    if (idx_z > NTotal)
                        break;
                }
                if (idx_z > NTotal)
                    break;
            }

            // --- 插值到接收/源深度（Sz, Rz）---
            for (int izz = 0; izz < NzTab; ++izz)
            {
                if (zTab(izz) > z_ssp(z_ssp.size() - 1))
                {
                    eigen.PhiMode(mode, izz) = 0.0;
                    dPhiMode(mode, izz) = 0.0;
                }
                else
                {
                    int idx = IZzTab(izz);
                    eigen.PhiMode(mode, izz) = std::complex<double>(Psi(idx)) +
                                               WTZ(izz) * std::complex<double>(dPsi(idx));
                    dPhiMode(mode, izz) = std::complex<double>(dPsidz(idx)) +
                                          WTZ(izz) * std::complex<double>(dPsidz(idx + 1) - dPsidz(idx));
                }
            }

            for (int isz = 0; isz < params.Pos.NSz; ++isz)
            {
                if (params.Pos.Sz(isz) > z_ssp(z_ssp.size() - 1))
                {
                    eigen.PsiS(mode, isz) = 0.0;
                    eigen.dPsidzS(mode, isz) = 0.0;
                }
                else
                {
                    if (NzTab == 1)
                    {
                        eigen.PsiS(mode, isz) = eigen.PhiMode(mode, 0);
                        eigen.dPsidzS(mode, isz) = dPhiMode(mode, 0);
                    }
                    else
                    {
                        int idx = ISzTab(isz);
                        eigen.PsiS(mode, isz) = eigen.PhiMode(mode, idx) +
                                                WTS(isz) * (eigen.PhiMode(mode, idx + 1) - eigen.PhiMode(mode, idx));
                        eigen.dPsidzS(mode, isz) = dPhiMode(mode, idx) +
                                                   WTS(isz) * (dPhiMode(mode, idx + 1) - dPhiMode(mode, idx));
                    }
                }
            }

            for (int irz = 0; irz < params.Pos.NRz; ++irz)
            {
                if (params.Pos.Rz(irz) > z_ssp(z_ssp.size() - 1))
                {
                    eigen.PsiR(mode, irz) = 0.0;
                    eigen.dPsidzR(mode, irz) = 0.0;
                }
                else
                {
                    if (NzTab == 1)
                    {
                        eigen.PsiR(mode, irz) = eigen.PhiMode(mode, 0);
                        eigen.dPsidzR(mode, irz) = dPhiMode(mode, 0);
                    }
                    else
                    {
                        int idx = IRzTab(irz);
                        eigen.PsiR(mode, irz) = eigen.PhiMode(mode, idx) +
                                                WTR(irz) * (eigen.PhiMode(mode, idx + 1) - eigen.PhiMode(mode, idx));
                        eigen.dPsidzR(mode, irz) = dPhiMode(mode, idx) +
                                                   WTR(irz) * (dPhiMode(mode, idx + 1) - dPhiMode(mode, idx));
                    }
                }
            }

            // --- 构建 phiZ 并转为 float（用于 .mod 文件）---
        }
        if (eigen.M > 0)
        {
            const double maxPsiS = eigen.PsiS.block(0, 0, eigen.M, eigen.PsiS.cols()).cwiseAbs().maxCoeff();
            const double maxPsiR = eigen.PsiR.block(0, 0, eigen.M, eigen.PsiR.cols()).cwiseAbs().maxCoeff();
            if (maxPsiS == 0.0 || maxPsiR == 0.0)
            {
                throw std::runtime_error("KRAKEN eigen solve produced zero source or receiver eigenfunctions.");
            }
        }
    }

    void Normalize(const size_t &iprof, const int &mode, const int &firstM, Eigen::VectorXd &Phi, int &ITP, int &NTotal1, double &x, TridMtx &trid, const OOK_parameters &params, EigenParams &eigen)
    {
        int iPower, modeCount = 0, j, j1, L, L1;
        double omega = (2 * pi * params.freqinfo.freq);
        double omega2 = SQ(omega);
        double x1, x2, DrhoDx = 0.0, DetaDx = 0.0, SqNorm = 0.0, RN, ScaleFactor, Slow = 0.0; // rhoMedium, rho_omega_h2,
        double rhoMedium, rho_omega_h2;
        std::complex<double> Perturbation_k = 0.0, Del, fTop1, gTop1, fTop2, gTop2, fBot1, gBot1, fBot2, gBot2;
        bool isTop = false, isComplex = false;
        auto ssp_prof = params.SSP.at(iprof);
        auto HSTop = ssp_prof.HSTop;
        auto HSBot = ssp_prof.HSBot;

        if (HSTop.BC == BC_Mode::MODE_A_Half_space)
        {
            std::complex<double> cp2 = SQ(HSTop.cp);
            Del = I1D * std::imag(HalfSpaceRadiationRoot(std::complex<double>(x - omega2 / cp2)));
            Perturbation_k -= Del * SQ(Phi(0)) / HSTop.rho;
            Slow += SQ(Phi(0)) / (2.0 * std::sqrt(x - std::real(omega2 / cp2))) / (HSTop.rho * std::real(cp2));
        }
        else if (HSTop.BC == BC_Mode::MODE_F_File || HSTop.BC == BC_Mode::MODE_P_Precomputed)
        {
            isTop = true;
            BCImpedance(iprof, x, isTop, fTop1, gTop1, iPower, isComplex, trid,
                        params, modeCount);
            isComplex = true;
            BCImpedance(iprof, x, isTop, fTop2, gTop2, iPower, isComplex, trid,
                        params, modeCount);
            Del = fTop2 / gTop2 - fTop1 / gTop1;
            Perturbation_k -= Del * SQ(Phi(0));
        }
        L = trid.Loc(ssp_prof.FirstAcoustic) - 1;
        j = 0;
        for (int im = ssp_prof.FirstAcoustic; im <= ssp_prof.LastAcoustic; im++)
        {
            // Compute contribution from the volume
            L += 1;
            rhoMedium = trid.rho(L);
            rho_omega_h2 = rhoMedium * omega2 * SQ(trid.h(im));
            // top interface
            SqNorm += 0.5 * trid.h(im) * SQ(Phi(j)) / rhoMedium;
            Slow += 0.5 * trid.h(im) * (trid.B1(L) + 2.0) * SQ(Phi(j)) / rho_omega_h2;
            Perturbation_k += 0.5 * trid.h(im) * I1D * trid.B1C(L) * SQ(Phi(j)) / rhoMedium;

            // medium
            L1 = L + 1;
            L += trid.N(im) - 1;
            j1 = j + 1;
            j += trid.N(im) - 1;

            SqNorm += trid.h(im) * Phi.segment(j1, j - j1 + 1).array().square().sum() / rhoMedium;
            Slow += trid.h(im) * ((trid.B1.segment(L1, L - L1 + 1).array() + 2.0) * Phi.segment(j1, j - j1 + 1).array().square()).sum() / rho_omega_h2;
            Perturbation_k += trid.h(im) * I1D * (trid.B1C.segment(L1, L - L1 + 1).array() * Phi.segment(j1, j - j1 + 1).array().square()).sum() / rhoMedium;

            // 底部边界：更新索引
            L += 1;
            j += 1;

            // 累加平方范数（SqNorm）
            SqNorm += 0.5 * trid.h(im) * (Phi(j) * Phi(j)) / rhoMedium;

            // 累加慢度项（Slow）
            Slow += 0.5 * trid.h(im) * (trid.B1(L) + 2.0) * SQ(Phi(j)) / rho_omega_h2;

            // 累加底部接口的扰动项（第二个Perturbation_k更新）
            Perturbation_k += 0.5 * trid.h(im) * I1D * trid.B1C(L) * SQ(Phi(j)) / rhoMedium;
        }

        if (HSBot.BC == BC_Mode::MODE_A_Half_space || HSBot.BC == BC_Mode::MODE_F_File || HSBot.BC == BC_Mode::MODE_P_Precomputed)
        {
            isTop = false;
            isComplex = false;
            BCImpedance(iprof, x, isTop, fBot1, gBot1, iPower, isComplex, trid,
                        params, modeCount);
            isComplex = true;
            BCImpedance(iprof, x, isTop, fBot2, gBot2, iPower, isComplex, trid,
                        params, modeCount);
            Del = fBot2 / gBot2 - fBot1 / gBot1;
            Perturbation_k -= Del * SQ(Phi(j));
            if (HSBot.BC == BC_Mode::MODE_A_Half_space)
            {
                std::complex<double> cp2 = SQ(HSBot.cp);
                Slow += SQ(Phi(j)) / (2.0 * std::sqrt(x - std::real(omega2 / cp2))) /
                        (HSBot.rho * std::real(cp2));
            }
        }

        // Compute derivative of top admitance
        x1 = 0.9999999 * x;
        x2 = 1.0000001 * x;
        isComplex = false;
        isTop = true;
        BCImpedance(iprof, x1, isTop, fTop1, gTop1, iPower, isComplex, trid,
                    params, modeCount);
        BCImpedance(iprof, x2, isTop, fTop2, gTop2, iPower, isComplex, trid,
                    params, modeCount);
        if (gTop1 != 0.0)
            DrhoDx = std::real((fTop2 / gTop2 - fTop1 / gTop1)) / (x2 - x1);

        // Compute derivative of bottom admitance
        isTop = false;
        BCImpedance(iprof, x1, isTop, fBot1, gBot1, iPower, isComplex, trid,
                    params, modeCount);
        BCImpedance(iprof, x2, isTop, fBot2, gBot2, iPower, isComplex, trid,
                    params, modeCount);
        if (gBot1 != 0.0)
            DetaDx = std::real((fBot2 / gBot2 - fBot1 / gBot1)) / (x2 - x1);

        // Scale the mode
        RN = SqNorm - DrhoDx * SQ(Phi(0)) + DetaDx * SQ(Phi(NTotal1 - 1));

        if (RN <= 0.0)
        {
            RN = -RN;
        }

        ScaleFactor = 1.0 / sqrt(RN);
        if (Phi(ITP) < 0.0)
            ScaleFactor = -ScaleFactor;

        Phi *= ScaleFactor;
        // std::cout << "Phi: " << Phi << std::endl;
        Slow = SQ(ScaleFactor) * Slow * omega / sqrt(x);
        Perturbation_k *= SQ(ScaleFactor);
        eigen.VG(mode) = 1 / Slow;
        ScatterLoss(iprof, mode, Perturbation_k, Phi, x, trid, params, eigen);
    }

    void ScatterLoss(const size_t &iprof, const int &mode, std::complex<double> &Perturbation_k, Eigen::VectorXd &Psi, double &x, TridMtx &trid, const OOK_parameters &params, EigenParams &eigen)
    {
        double omega = 2 * pi * params.freqinfo.freq;
        double omega2 = SQ(omega);
        double rho1, rho2, rhoInside, h2, sigma;
        std::complex<double> eta1Sq, eta2Sq, U, PhiC; // kx = sqrt(x),
        auto ssp_prof = params.SSP.at(iprof);
        auto HSTop = ssp_prof.HSTop;
        auto HSBot = ssp_prof.HSBot;

        int j = 0, L = trid.Loc(ssp_prof.FirstAcoustic);
        for (int im = ssp_prof.FirstAcoustic; im <= ssp_prof.LastAcoustic; im++) // Loop over media
        {
            // Calculate rho1, eta1Sq, Phi, U
            if (im == ssp_prof.FirstAcoustic) // Top properties
            {
                switch (HSTop.BC)
                {
                case BC_Mode::MODE_A_Half_space:
                    rho1 = HSTop.rho;
                    eta1Sq = x - omega2 / SQ(HSTop.cp);
                    U = sqrt(eta1Sq) * Psi(0) / rho1;
                    break;
                case BC_Mode::MODE_V_Vacuum:
                    rho1 = 1e-9;
                    eta1Sq = 1.0;
                    rhoInside = trid.rho(trid.Loc(params.SSP[iprof].FirstAcoustic));
                    U = Psi(1) / trid.h(params.SSP[iprof].FirstAcoustic) / rhoInside;
                    break;
                case BC_Mode::MODE_R_Rigid:
                    rho1 = 1e9;
                    eta1Sq = 1.0;
                    U = 0.0;
                    break;
                default:
                    rho1 = 0.0;
                    eta1Sq = 0.0;
                    U = 0.0;
                    break;
                }
                sigma = ssp_prof.sigma[im];
            }
            else
            {
                h2 = SQ(trid.h(im));
                j = j + trid.N(im);
                L = trid.Loc(im) + trid.N(im);
                rho1 = trid.rho(L);
                eta1Sq = (2.0 + trid.B1(L)) / h2 - x;
                U = (-Psi(j - 1) - 0.5 * (trid.B1(L) - h2 * x) * Psi(j)) / (trid.h(im) * rho1);
                sigma = ssp_prof.sigma[im];
            }

            // Calculate rho2, eta2
            if (im == params.SSP[iprof].LastAcoustic) // Bottom properties
            {
                switch (HSBot.BC)
                {
                case BC_Mode::MODE_A_Half_space: // Acousto-elastic
                    rho2 = HSBot.rho;
                    eta2Sq = omega2 / SQ(HSBot.cp) - x;
                    break;
                case BC_Mode::MODE_V_Vacuum: // Vacuum
                    rho2 = 1e-9;
                    eta2Sq = 1.0;
                    break;
                case BC_Mode::MODE_R_Rigid: // Rigid
                    rho2 = 1e9;
                    eta2Sq = 1.0;
                    break;
                default: // Tabulated
                    rho2 = 0.0;
                    eta2Sq = 0.0;
                    break;
                }
                sigma = ssp_prof.sigma[im];
            }
            else
            {
                rho2 = trid.rho(L + 1);
                eta2Sq = (2.0 + trid.B1(L + 1)) / SQ(trid.h(im + 1)) - x;
            }
            PhiC = Psi(j);
            std::complex<double> kup = KupIng(sigma, eta1Sq, rho1, eta2Sq, rho2, PhiC, U);
            // // 打印kup
            // cout << "kup: \n"
            //      << kup << endl;
            Perturbation_k += kup;
        }
        // 打印Perturbation_k
        // cout << "Perturbation_k: \n"
        //      << Perturbation_k << endl;
        eigen.k(mode) = Perturbation_k;
    }
}
