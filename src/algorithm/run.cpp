#include "run.h"

void run()
{
    parameters params;
    set_pekeris(params);
    // set_Munk(params);

    // 检查竖直网络步长是否小于波长/20
    for (int i = 0; i < params.NMedia; ++i)
    {
        double h = params.SSP[i].depth / params.SSP[i].N;
        double lambda_1_20 = params.SSP[i].alphaR[params.SSP[i].NPts - 1] / params.freqinfo->freq / 20.0; // 最后一个声速计算波长
        int Nneeded = int((params.SSP[i].depth) / lambda_1_20);
        Nneeded = std::max(Nneeded, 10); // require a minimum of 10 points				要求每一层媒质至少有10个点

        if (h > lambda_1_20)
        {
            // 打印警告信息（中文）
            cout << "警告：KRAKEN 垂直网格步长太大，已经将网格数量: " << params.SSP[i].N << " 调整为: " << Nneeded << endl;
            params.SSP[i].N = Nneeded;
        }
        else if (params.mesh.N(i) == 0) // 网格数为0时，将网格数调整为Nneeded
        {
            params.SSP[i].N = Nneeded;
        }
    }
    // 初始化cp cs
    UpdateSSPLoss(params.freqinfo->freq, params.freqinfo->freq, params.NMedia,
                  params.SSPType, params.AttenUnit, params.SSP);
    UpdateHSLoss(params.freqinfo->freq, params.freqinfo->freq, params.NMedia,
                 params.AttenUnit, params.HSTop, params.HSBot);

    // TODO 计算本征值和本征函数
    KrakenMatrix kramtrx;
    EigenParams eigen;
    EigenFunction eigenfun;
    double error;
    params.mesh.hV.resize(params.mesh.NSets);
    for (int iset = 0; iset < params.mesh.NSets; iset++)
    {
        int ntimes = params.mesh.NV[iset];
        Initialize(iset, params, kramtrx, ntimes);
        SolveEp(iset, params.mesh.NSets, eigen, eigenfun, kramtrx, params, error);
        if (error * params.Rmax < 1.0)
        {
            break;
        }
        else
        {
            if (iset == params.mesh.NSets - 1)
                cout << "Warning in KRAKEN : Too many meshes needed: check convergence" << endl;
        }
    }

    int M = 0;
    VectorXd Ex1 = eigen.Extrap.segment(0, eigen.M);
    // cout << "Ex1: \n"
    //      << Ex1 << endl;
    while (M < eigen.M && Ex1(M) > SQ(2 * pi * params.freqinfo->freq / params.Chigh))
    {
        M++;
    }
    eigen.M = M;
    // 使用tempk存储本征值
    VectorXcd tempk = eigen.k.segment(0, M);
    eigen.k.resize(M);
    eigen.k = tempk;
    // 使用tempphi存储本征函数
    MatrixXcd tempphi = eigenfun.phi.block(0, 0, eigen.M, eigenfun.phi.cols());
    eigenfun.phi.resize(M, eigenfun.phi.cols());
    eigenfun.phi = tempphi;
    // 使用tempphiR存储本征函数的实部
    MatrixXcd tempphiR = eigenfun.phiR.block(0, 0, eigen.M, eigenfun.phiR.cols());
    eigenfun.phiR.resize(M, eigenfun.phiR.cols());
    eigenfun.phiR = tempphiR;
    // 使用tempphiS存储本征函数的虚部
    MatrixXcd tempphiS = eigenfun.phiS.block(0, 0, eigen.M, eigenfun.phiS.cols());
    eigenfun.phiS.resize(M, eigenfun.phiS.cols());
    eigenfun.phiS = tempphiS;
    // 使用tempdphidz存储本征函数的导数
    MatrixXcd tempdphidz = eigenfun.dphidz.block(0, 0, eigen.M, eigenfun.dphidz.cols());
    eigenfun.dphidz.resize(M, eigenfun.dphidz.cols());
    eigenfun.dphidz = tempdphidz;
    // 使用tempdphidzR存储本征函数的导数的实部
    MatrixXcd tempdphidzR = eigenfun.dphidzR.block(0, 0, eigen.M, eigenfun.dphidzR.cols());
    eigenfun.dphidzR.resize(M, eigenfun.dphidzR.cols());
    eigenfun.dphidzR = tempdphidzR;
    // 使用tempdphidzS存储本征函数的导数的虚部
    MatrixXcd tempdphidzS = eigenfun.dphidzS.block(0, 0, eigen.M, eigenfun.dphidzS.cols());
    eigenfun.dphidzS.resize(M, eigenfun.dphidzS.cols());
    eigenfun.dphidzS = tempdphidzS;

    // // 打印eigen.k
    // cout << "eigen.k: \n"
    //      << eigen.k << endl;
    for (int i = 0; i < M; i++)
    {
        eigen.k(i) = sqrt(eigen.Extrap(i) + eigen.k(i));
    }
    // // 打印eigen.k
    // cout << "eigen.k: \n"
    //      << eigen.k << endl;

    // 写入.mod
    // 将 Fortran 的 WRITE 逻辑转换为 C++ 文件写入
    // 假设 MODFile 已以二进制方式打开，且已定位到 IRecProfile 对应位置
    // 此处用 std::ofstream 模拟，实际工程中请确保文件已正确打开并保持同步

    // 写入模态数 M
    params.MODFile.seekp(eigen.IRecProfile * 4 * eigen.LRecordLength, std::ios::beg);
    params.MODFile.write(reinterpret_cast<const char *>(&M), sizeof(int));

    // 写入复本征值 k
    int IFirst = 0; // C++ 从 0 开始
    for (int IREC = 0; IREC < (2 * M - 1) / eigen.LRecordLength + 1; ++IREC)
    {
        int ILast = std::min(M, IFirst + eigen.LRecordLength / 2) - 1;
        int segLen = ILast - IFirst + 1;

        // 定位到对应记录
        params.MODFile.seekp((eigen.IRecProfile + 2 + M + IREC) * 4 * eigen.LRecordLength, std::ios::beg);

        // 将 eigen.k 中 IFirst 到 ILast 的复数写出
        // 转换为complex float输出

        for (int i = 0; i < segLen; i++)
        {
            complex<float> kf;
            kf = std::complex<float>(eigen.k(IFirst + i).real(), eigen.k(IFirst + i).imag());
            params.MODFile.write(reinterpret_cast<const char *>(&kf),
                                 sizeof(std::complex<float>));
        }

        IFirst = ILast + 1;
    }

    // 更新下一段起始记录号
    eigen.IRecProfile += 3 + M + (2 * M - 1) / eigen.LRecordLength;

    params.MODFile.close();

    size_t N = (size_t)params.Pos->NSz * (size_t)params.Pos->NRz_per_range * (size_t)params.Pos->NRr;
    std::complex<float> *u_AllSources;
    u_AllSources = new std::complex<float>[N];
    std::complex<float> *uAllSources_vr;
    uAllSources_vr = new std::complex<float>[N];
    std::complex<float> *uAllSources_vz;
    uAllSources_vz = new std::complex<float>[N];
    for (int isz = 0; isz < params.Pos->NSz; isz++)
    {
        Evaluate(eigenfun, eigen, params, isz, u_AllSources, uAllSources_vr, uAllSources_vz);
    }
    string filename = "test_pressure";
    string filename_vr = "test_vr";
    string filename_vz = "test_vz";
    export_shd(filename, params, u_AllSources);
    export_shd(filename_vr, params, uAllSources_vr);
    export_shd(filename_vz, params, uAllSources_vz);
}
