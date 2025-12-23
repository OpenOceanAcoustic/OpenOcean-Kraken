
#include "paramsBase.h"

class input_SSP : public paramsBase
{
public:
    input_SSP() {}
    virtual ~input_SSP() {}

    // 初始化函数
    virtual void Init(parameters &params) const override
    {   
        params.SSP = nullptr;// 先初始化为空
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

        for (int iprof = 0; iprof < params.NProf; iprof++)
        {
            auto &ssp = params.SSP[iprof];
            // 先初始化内存
            int vec_size = ssp.NPts.sum();
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
    void set_SSP(parameters &params, const std::vector<SSP_1D>& sspVec) const
    {   
        
        // 1. 先清理旧内存（防止泄漏）
        params.SSP.reset(); //清理智能指针
        params.NProf = 0;

        size_t N = sspVec.size();
        if (N == 0) return;

        params.SSP = std::make_unique<SSPStructure[]>(N);//只能指针
        params.NProf = N;
        for (size_t i = 0; i < N; ++i) {
            params.SSP[i] = convert_to_SSPStructure(sspVec[i]); // 调用拷贝赋值运算符
        }
    }

    void set_AttenUnit(parameters &params, Atten_Mode unit) const
    {
        params.AttenUnit = unit;
    }

    SSPStructure convert_to_SSPStructure(const SSP_1D& ssp1d) const {
        SSPStructure dst;
        const auto& layers = ssp1d.layers;
        int nmedia = static_cast<int>(layers.size());

        if (nmedia == 0) {
            dst.NMedia = 0;
            dst.FirstAcoustic = -1;
            dst.LastAcoustic = -1;
            return dst;
        }

        // 1. 基本参数（从 ssp1d 继承）
        dst.NMedia = nmedia;
        dst.SSPType = ssp1d.SSPType;

        // 假设所有层都是声学层（按需调整）
        dst.FirstAcoustic = 0;
        dst.LastAcoustic = nmedia - 1;

        // 2. 分配层元数据
        dst.NPts.resize(nmedia);
        dst.NMesh.resize(nmedia);
        dst.beta.resize(nmedia);
        dst.ft.resize(nmedia);
        dst.sigma.resize(nmedia);
        dst.offset.resize(nmedia);
        dst.Material.resize(nmedia);

        // 注意：dst.depth 可能不需要！除非你明确要用
        // 如果必须填，可填每层最大 z（底部深度）或忽略
        dst.depth.resize(nmedia); // 暂时保留，按需使用

        // 3. 计算 offset 和总点数
        int total_pts = 0;
        for (int i = 0; i < nmedia; ++i) {
            dst.offset[i] = total_pts;
            const auto& layer = layers[i];
            dst.NPts[i] = layer.npts;
            dst.NMesh[i] = layer.nmesh;
            dst.beta[i] = layer.beta;
            dst.ft[i] = layer.ft;
            dst.sigma[i] = layer.sigma;
            dst.Material[i] = layer.Material;
            // 可选：填 depth 为该层最大深度（假设 z 单调增）
            if (layer.npts > 0) {
                dst.depth[i] = layer.z.maxCoeff(); // 或 layer.z(layer.npts - 1)
            } else {
                dst.depth[i] = (i == 0) ? 0.0 : dst.depth[i - 1];
            }

            total_pts += layer.npts;
        }

        // 4. 利用 flatten() 避免重复代码！
        auto flat = ssp1d.flatten();

        // 直接赋值（避免手动 segment 循环）
        dst.z = std::move(flat.z);
        dst.rho = std::move(flat.rho);
        dst.alphaR = std::move(flat.alphaR);
        dst.alphaI = std::move(flat.alphaI);
        dst.betaR = std::move(flat.betaR);
        dst.betaI = std::move(flat.betaI);

        // 6. 清空插值系数（由后续函数填充）
        dst.row_size = 4;
        dst.cspline.resize(0, 0);
        dst.cpCoef.resize(0, 0);
        dst.csCoef.resize(0, 0);
        dst.rhoCoef.resize(0, 0);
        dst.cpSpline.resize(0, 0);
        dst.csSpline.resize(0, 0);
        dst.rhoSpline.resize(0, 0);
        dst.cCoef.resize(0, 0);
        dst.csWork.resize(0, 0);

        return dst;
    }




    // 打印
    virtual void Echo(parameters &params) const override
    {
    }

private:
    // 设置默认值的私有方法
    void setDefaultValues(parameters &params) const
    {   


        SSP_1D ssp1d;
        SSPLayer layer;


        layer.alphaR = VectorXd::Ones(2)*1500;
        layer.alphaI = VectorXd::Zero(2); // 纵波速度虚部（或衰减）
        layer.betaR = VectorXd::Zero(2);
        layer.betaI = VectorXd::Zero(2);
        layer.rho = VectorXd::Ones(2);
        layer.z = Vector2d(0, 200.0);
        layer.nmesh = 0;
        layer.sigma = 0;
        layer.npts = 2;
        layer.Material = Media_Mode::MODE_A_Acoustic;
        layer.beta = 0.0;
        layer.ft = 0.0;
        layer.npts = 2;

        ssp1d.add(layer);//添加层

        SSPStructure ssp = convert_to_SSPStructure(ssp1d);

        params.SSP = std::make_unique<SSPStructure[]>(1);//只能指针
        params.SSP[0] = ssp;
        params.NProf = 1;
 
  
    }
};
