
#include "paramsBase.h"

class input_SSP : public paramsBase
{
public:
    input_SSP() {}
    virtual ~input_SSP() {}

    // 初始化函数
    virtual void Init(parameters &params) const override
    {
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
        VectorXi NMeshMaxVec(params.NProf);
        VectorXi NMediaVec(params.NProf);
        NMeshMaxVec.setZero();

        for (size_t iprof = 0; iprof < params.NProf; iprof++)
        {
            auto &ssp = params.SSP[iprof];
            // 先初始化内存
            size_t vec_size = ssp.NPts.sum();
            ssp.cp.resize(vec_size);
            ssp.cs.resize(vec_size);
            if (ssp.SSPType == SSP_Mode::MODE_P_cPCHIP)
            {
                ssp.cpCoef.resize(ssp.row_size, vec_size);
                ssp.csCoef.resize(ssp.row_size, vec_size);
                ssp.rhoCoef.resize(ssp.row_size, vec_size);
            }
            if (ssp.SSPType == SSP_Mode::MODE_S_cCubic)
            {
                ssp.cpSpline.resize(ssp.row_size, vec_size);
                ssp.csSpline.resize(ssp.row_size, vec_size);
                ssp.rhoSpline.resize(ssp.row_size, vec_size);
            }

            UpdateSSPLoss(freq, freq, ssp.NMedia, ssp.SSPType, params.AttenUnit, ssp);
            ssp.offset[0] = 0;
            ssp.FirstAcoustic = -1;
            size_t ilay = 0;
            size_t NMeshAll = 0;
            NMediaVec[iprof] = ssp.NMedia;
            ssp.Material.resize(ssp.NMedia);
            for (int imedia = 0; imedia < ssp.NMedia; imedia++)
            {
                if (imedia > 0)
                    ssp.offset[imedia] = ssp.offset[imedia - 1] + ssp.NPts[imedia - 1];
                if (std::real(ssp.cs[0]) == 0.0)
                { // 声学介质情况
                    ssp.Material[imedia] = Media_Mode::MODE_A_Acoustic;
                    if (ssp.FirstAcoustic == -1)
                    {
                        ssp.FirstAcoustic = imedia;
                    }
                }
                ssp.LastAcoustic = imedia;
                ilay += ssp.NPts[imedia];
                double h = ssp.depth[imedia] / ssp.NMesh[imedia];
                double lambda_1_20 = ssp.alphaR[ilay - 1] / params.freqinfo->freq / 20.0; // 最后一个声速计算波长
                int Nneeded = int((ssp.depth[imedia]) / lambda_1_20);
                Nneeded = std::max(Nneeded, 10); // require a minimum of 10 points				要求每一层媒质至少有10个点

                if (ssp.NMesh[imedia] == 0) // 网格数为0时，将网格数调整为Nneeded
                {
                    // 打印信息（中文）
                    cout << "网格数为0，按照波长1/20计算，已经将网格数量设定为: " << Nneeded << endl;
                    ssp.NMesh[imedia] = Nneeded;
                }
                else if (h > lambda_1_20)
                {
                    // 打印警告信息（中文）
                    cout << "警告：KRAKEN 垂直网格步长太大，已经将网格数量: " << ssp.NMesh[imedia] << " 调整为: " << Nneeded << endl;
                    ssp.NMesh[imedia] = Nneeded;
                }
                NMeshAll += ssp.NMesh[imedia];
            }
            // 每一层都需要+1
            NMeshAll += ssp.NMedia;
            NMeshMaxVec[iprof] = NMeshAll * params.mesh.NV[1]; // iset=1时最大的网格数
        }
        params.NMeshMax = NMeshMaxVec.maxCoeff(); // 所有NProf中最大的网格数
        params.NMediaMax = NMediaVec.maxCoeff(); // 所有NProf中最大的媒质数
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
        params.SSP[0].beta[0] = 0.0;
        params.SSP[0].ft[0] = 0.0;
        params.SSP[0].NMesh[0] = 0;
        params.SSP[0].offset[0] = 0;
        params.SSP[0].depth[0] = 200;
        params.SSP[0].sigma[0] = 0;
    }
};
