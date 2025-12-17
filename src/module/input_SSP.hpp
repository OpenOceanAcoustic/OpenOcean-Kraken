
#include "paramsBase.h"

class input_SSP : public paramsBase
{
public:
    input_SSP() {}
    virtual ~input_SSP() {}

    // 初始化函数
    virtual void Init(parameters &params) const override
    {
        params.SSP = new SSPStructure[params.NProf];
    }

    // 设置默认值
    virtual void Default(parameters &params) const override
    {
        setDefaultValues(params);
    }

    // 运行前的预处理
    virtual void Preprocess(parameters &params) const override
    {
        double freq = params.freqinfo->freq;
        for (size_t iprof = 0; iprof < params.NProf; iprof++)
        {
            auto &ssp = params.SSP[iprof];
            ssp.offset[0] = 0;
            size_t ilay = 0;
            size_t NMeshAll = 0;
            for (int imedia = 0; imedia < ssp.NMedia; imedia++)
            {
                if (imedia > 0)
                    ssp.offset[imedia] = ssp.offset[imedia - 1] + ssp.NPts[imedia - 1];
                ilay += params.SSP[iprof].NPts[imedia];
                double h = params.SSP[iprof].depth[imedia] / params.SSP[iprof].NMesh[imedia];
                double lambda_1_20 = params.SSP[iprof].alphaR[ilay - 1] / params.freqinfo->freq / 20.0; // 最后一个声速计算波长
                int Nneeded = int((params.SSP[iprof].depth[imedia]) / lambda_1_20);
                Nneeded = std::max(Nneeded, 10); // require a minimum of 10 points				要求每一层媒质至少有10个点

                if (h > lambda_1_20)
                {
                    // 打印警告信息（中文）
                    cout << "警告：KRAKEN 垂直网格步长太大，已经将网格数量: " << params.SSP[iprof].NMesh[imedia] << " 调整为: " << Nneeded << endl;
                    params.SSP[iprof].NMesh[imedia] = Nneeded;
                }
                else if (params.mesh.N(imedia - 1) == 0) // 网格数为0时，将网格数调整为Nneeded
                {
                    params.SSP[iprof].NMesh[imedia] = Nneeded;
                }
                NMeshAll += params.SSP[iprof].NMesh[imedia];
            }

            params.SSP[iprof].FirstAcoustic = -1;
            params.mesh[iprof].Loc.resize(params.SSP[iprof].NMedia + 1);
            params.mesh[iprof].Loc[0] = 0; // C++使用0-based索引

            // 计算总网格点数
            params.mesh[iprof].N.resize(params.SSP[iprof].NMedia);
            params.mesh[iprof].h.resize(params.SSP[iprof].NMedia);

            // 每一层都需要+1
            NMeshAll += ssp.NMedia;
            ssp.NMeshMax = NMeshAll * params.mesh.NV[params.mesh.NSets - 1]; // 最大的网格数
            // 最大内存初始化，避免每次的维度不一致
            ssp.cp_int.resize(ssp.NMeshMax);
            ssp.cs_int.resize(ssp.NMeshMax);
            ssp.rho_int.resize(ssp.NMeshMax);

            // 将计算矩阵也初始化
            ssp.B1.resize(ssp.NMeshMax);
            ssp.B1C.resize(ssp.NMeshMax);
            ssp.B2.resize(ssp.NMeshMax);
            ssp.B3.resize(ssp.NMeshMax);
            ssp.B4.resize(ssp.NMeshMax);
            ssp.rhoparam.resize(ssp.NMeshMax);
            if (ssp.SSPType == SSP_Mode::MODE_P_cPCHIP)
            {
                ssp.cpCoef.resize(4, ssp.NMeshMax);
                ssp.csCoef.resize(4, ssp.NMeshMax);
                ssp.rhoCoef.resize(4, ssp.NMeshMax);
            }
            else if (ssp.SSPType == SSP_Mode::MODE_S_cCubic)
            {
                ssp.cpSpline.resize(4, ssp.NMeshMax);
                ssp.csSpline.resize(4, ssp.NMeshMax);
                ssp.rhoSpline.resize(4, ssp.NMeshMax);
            }
            // 计算SSP的参数
            UpdateSSPLoss(freq, freq, ssp.NMedia,
                          ssp.SSPType, params.AttenUnit, ssp);
        }
    }

    // 设置SSP的方法
    void set_SSP(parameters &params, SSP_1D *sspInput, size_t NProf) const
    {
        params.SSP = new SSPStructure[NProf];
        for (size_t iprof = 0; iprof < NProf; iprof++)
        {
            auto &ssp = params.SSP[iprof];
            ssp.NMedia = sspInput[iprof].NMedia;
            ssp.NPts = sspInput[iprof].NPts;
            ssp.NMesh = sspInput[iprof].NMesh;
            ssp.SSPType = sspInput[iprof].SSPType;
            ssp.beta = sspInput[iprof].beta;
            ssp.ft = sspInput[iprof].ft;
            ssp.sigma = sspInput[iprof].sigma;
            ssp.z = sspInput[iprof].z;
            ssp.rho = sspInput[iprof].rho;
            ssp.alphaR = sspInput[iprof].alphaR;
            ssp.alphaI = sspInput[iprof].alphaI;
            ssp.betaR = sspInput[iprof].betaR;
            ssp.betaI = sspInput[iprof].betaI;
        }
    }

    void set_AttenUnit(parameters &params, Atten_Mode unit) const
    {
        params.AttenUnit = unit;
    }

    SSP_1D get_SSP(const parameters &params, size_t iprof) const
    {
        auto &ssp = params.SSP[iprof];
        SSP_1D sspOutput;
        sspOutput.NMedia = ssp.NMedia;
        sspOutput.NPts = ssp.NPts;
        sspOutput.NMesh = ssp.NMesh;
        sspOutput.SSPType = ssp.SSPType;
        sspOutput.beta = ssp.beta;
        sspOutput.ft = ssp.ft;
        sspOutput.sigma = ssp.sigma;
        sspOutput.z = ssp.z;
        sspOutput.rho = ssp.rho;
        sspOutput.alphaR = ssp.alphaR;
        sspOutput.alphaI = ssp.alphaI;
        sspOutput.betaR = ssp.betaR;
        sspOutput.betaI = ssp.betaI;
        return sspOutput;
    }

    // 打印
    virtual void Echo(parameters &params) const override
    {
    }

private:
    // 设置默认值的私有方法
    void setDefaultValues(parameters &params) const
    {
        params.SSP[0].NMedia = 1;
        // 声速剖面类型
        params.SSP[0].SSPType = SSP_Mode::MODE_C_cLinear;

        // 声速剖面
        params.SSP[0].NPts.resize(params.SSP[0].NMedia);
        params.SSP[0].beta.resize(params.SSP[0].NMedia);
        params.SSP[0].ft.resize(params.SSP[0].NMedia);
        params.SSP[0].NMesh.resize(params.SSP[0].NMedia);
        params.SSP[0].depth.resize(params.SSP[0].NMedia);
        params.SSP[0].sigma.resize(params.SSP[0].NMedia);
        params.SSP[0].offset.resize(params.SSP[0].NMedia);

        params.SSP[0].NPts[0] = 2;
        params.SSP[0].alphaR = Vector2d(1500.0, 1500.0);
        params.SSP[0].alphaI = Vector2d(0.0, 0.0);
        params.SSP[0].betaR = Vector2d(0.0, 0.0);
        params.SSP[0].betaI = Vector2d(0.0, 0.0);
        params.SSP[0].rho = Vector2d(1.0, 1.0);
        params.SSP[0].z = Vector2d(0, 200.0);
        params.SSP[0].cp.resize(params.SSP[0].NPts[0]);
        params.SSP[0].cs.resize(params.SSP[0].NPts[0]);
        params.SSP[0].beta[0] = 0.0;
        params.SSP[0].ft[0] = 0.0;
        params.SSP[0].NMesh[0] = 0;
        params.SSP[0].offset[0] = 0;
        params.SSP[0].depth[0] = 200;
        params.SSP[0].sigma[0] = 0;
    }
};
