#include "paramsBase.h"
#include "AttenMod.h"
#include <stdexcept>
#include <utility>
namespace OpenOceanKraken
{
    class input_SSP : public paramsBase
    {
    public:
        input_SSP() {}
        virtual ~input_SSP() {}

        // 初始化函数
        virtual void Init(OOK_parameters &params) const override
        {
        }

        // 设置默认值
        virtual void Default(OOK_parameters &params) const override
        {
            setMunk(params);
        }

        // 运行前的预处理
        virtual void Preprocess(OOK_parameters &params) const override
        {
            int NSSP = params.SSP.size();
            double freq = params.freqinfo.freq;
            Eigen::VectorXi NMeshMaxVec(NSSP);
            Eigen::VectorXi NMediaVec(NSSP);
            NMeshMaxVec.setZero();

            for (int iprof = 0; iprof < NSSP; iprof++)
            {
                auto &ssp = params.SSP.at(iprof);
                // === 1. 分配 cp, cs 等扁平数组内存 ===
                int vec_size = ssp.NPts.sum();
                ssp.cp.resize(vec_size);
                ssp.cs.resize(vec_size);
                // === 2. 根据 SSPType 分配插值系数内存 ===
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
                // === 3. 计算介质损耗（声吸收）===
                ssp::UpdateSSPLoss(freq, freq, ssp.NMedia, ssp.SSPType, params.AttenUnit, ssp); // 计算ssp损失
                ssp::UpdateHSLoss(freq, freq, params.AttenUnit, ssp.HSTop, ssp.HSBot);          // 计算边界损失

                // === 4. 初始化介质类型和声学层范围 ===
                ssp.Material.resize(ssp.NMedia);
                ssp.FirstAcoustic = -1;
                ssp.LastAcoustic = -1;

                size_t NMeshAll = 0;
                NMediaVec(iprof) = ssp.NMedia;

                for (int imedia = 0; imedia < ssp.NMedia; imedia++)
                {
                    // 获取当前层在扁平数组中的起始索引
                    size_t layer_start = ssp.get_media_start(imedia);
                    size_t layer_end = ssp.get_media_end(imedia);
                    // 判断是否为声学层：检查该层第一个点的 cs 是否为 0
                    if (std::real(ssp.cs[layer_start]) == 0.0)
                    { // 声学介质情况
                        ssp.Material[imedia] = Media_Mode::MODE_A_Acoustic;
                        if (ssp.FirstAcoustic == -1)
                        {
                            ssp.FirstAcoustic = imedia;
                        }
                        ssp.LastAcoustic = imedia;
                    }
                    else{
                        ssp.Material[imedia] = Media_Mode::MODE_E_Elastic;
                    }
                    

                    double lambda_1_20 = ssp.alphaR[layer_end - 1] / params.freqinfo.freq / 20.0; // 最后一个声速计算波长
                    if (lambda_1_20 <= 0.0)
                    {
                        throw std::runtime_error("Invalid SSP mesh estimate: frequency and sound speed must be positive.");
                    }
                    int Nneeded = int((ssp.depth[imedia]) / lambda_1_20);
                    if (ssp.NMesh[imedia] < 0)
                    {
                        throw std::runtime_error("Invalid SSP mesh count: NMesh must be non-negative.");
                    }
                    Nneeded = std::max(Nneeded, 10); // require a minimum of 10 points				要求每一层媒质至少有10个点

                    if (ssp.NMesh[imedia] == 0) // 网格数为0时，将网格数调整为Nneeded
                    {
                        // 打印信息（中文）
                        std::cout << "网格数为0，按照波长1/20计算，已经将网格数量设定为: " << Nneeded << std::endl;
                        ssp.NMesh[imedia] = Nneeded;
                    }
                    if (ssp.NMesh[imedia] < Nneeded / 2)
                    {
                        // 打印警告信息（中文）
                        std::cout << "警告：KRAKEN 垂直网格步长太大，已经将网格数量: " << ssp.NMesh[imedia] << " 调整为: " << Nneeded << std::endl;
                        throw std::runtime_error("SSP mesh count is too coarse.");
                    }
                    NMeshAll += ssp.NMesh[imedia];
                }
                // 每一层都需要+1
                NMeshAll += ssp.NMedia;
                int maxNV = params.mesh.NV[0];
                for (int iset = 1; iset < params.mesh.NSets && iset < 5; ++iset)
                {
                    maxNV = std::max(maxNV, params.mesh.NV[iset]);
                }
                NMeshMaxVec(iprof) = NMeshAll * maxNV;
            }
            params.NMeshMax = NMeshMaxVec.maxCoeff(); // 所有NProf中最大的网格数
            params.NMediaMax = NMediaVec.maxCoeff();  // 所有NProf中最大的媒质数
            std::cout<<"OpenOcean-Kraken: input SSP processed"<<std::endl;
        }

        // 设置SSP的方法
        void set_SSP(OOK_parameters &params, const std::vector<ssp::Range_Independent_Area> &sspVec) const
        {

            // 1. 先清理旧内存（防止泄漏）
            params.SSP.clear(); // 清理智能指针
            params.sspInput = sspVec; // 赋值一份
            size_t N = sspVec.size();
            if (N == 0)
                return;

            params.SSP.resize(N); //

            for (size_t i = 0; i < N; ++i)
            {
                params.SSP[i] = convert_to_SSPStructure(sspVec[i]); // 调用拷贝赋值运算符
            }
        }

        void set_AttenUnit(OOK_parameters &params, Atten_Mode unit) const
        {
            validateAttenuationMode(unit);
            params.AttenUnit = std::move(unit);
        }

        ssp::SSPStructure convert_to_SSPStructure(const ssp::Range_Independent_Area &Area_1D) const
        {
            ssp::SSPStructure dst;
            const auto &layers = Area_1D.layers;
            int nmedia = static_cast<int>(layers.size());

            if (nmedia == 0)
            {
                dst.NMedia = 0;
                dst.FirstAcoustic = -1;
                dst.LastAcoustic = -1;
                return dst;
            }

            // 1. 基本参数（从 ssp1d 继承）
            dst.NMedia = nmedia;
            dst.SSPType = Area_1D.SSPType;

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
            dst.interp_offset.resize(nmedia);
            // 注意：dst.depth 可能不需要！除非你明确要用
            // 如果必须填，可填每层最大 z（底部深度）或忽略
            dst.depth.resize(nmedia); // 暂时保留，按需使用

            // 3. 计算 offset 和总点数
            int total_pts = 0;
            int interp_off = 0;
            for (int i = 0; i < nmedia; ++i)
            {
                dst.offset[i] = total_pts;
                dst.interp_offset[i] = interp_off;
                const auto &layer = layers[i];
                dst.NPts[i] = layer.npts;
                dst.NMesh[i] = layer.nmesh;
                dst.beta[i] = layer.beta;
                dst.ft[i] = layer.ft;
                dst.sigma[i] = layer.sigma;
                dst.Material[i] = layer.Material;
                // 可选：填 depth 为该层最大深度（假设 z 单调增）
                if (layer.npts > 0)
                {
                    dst.depth[i] = layer.z.maxCoeff() - layer.z.minCoeff();
                }
                else
                {
                    dst.depth[i] = (i == 0) ? 0.0 : dst.depth[i - 1];
                }

                total_pts += layer.npts;
                interp_off += dst.NMesh[i] + 1;
            }

            // 4. 利用 flatten() 避免重复代码！
            auto flat = Area_1D.flatten();

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

            // 7. 底部与顶部半空间属性
            dst.HSTop = Area_1D.HSTop;
            dst.HSBot = Area_1D.HSBot;
            return dst;
        }

        // 打印
        virtual void Echo(OOK_parameters &params) const override
        {
        }

    private:
        // 设置默认值的私有方法
        void setMunk(OOK_parameters &params) const
        {

            ssp::Range_Independent_Area ssp1d; // 1D 范围独立层结构

            Eigen::VectorXd alphaR(27);
            Eigen::VectorXd z(27);
            Eigen::VectorXd alphaI = Eigen::VectorXd::Zero(27);
            Eigen::VectorXd betaR = Eigen::VectorXd::Zero(27);
            Eigen::VectorXd betaI = Eigen::VectorXd::Zero(27);
            Eigen::VectorXd rho = Eigen::VectorXd::Ones(27);
            z << 0, 200, 250, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600, 3800, 4000, 4200, 4400, 4600, 4800, 5000;
            alphaR << 1548.52, 1530.29, 1526.69, 1517.78, 1509.49, 1504.3, 1501.38, 1500.14, 1500.12, 1501.02, 1502.57, 1504.62, 1507.02, 1509.69, 1512.55, 1515.56, 1518.67, 1521.85, 1525.1, 1528.38, 1531.7, 1535.04, 1538.39, 1541.76, 1545.14, 1548.52, 1551.91;

            ssp::SSPLayer layer{27, 0, 0, 0, 0.0, z, alphaR, alphaI, betaR, betaI, rho, Media_Mode::MODE_A_Acoustic};

            // 设置ssp
            ssp1d.addLayer(layer); // 添加层
            ssp1d.SSPType = SSP_Mode::MODE_C_cLinear;
            // 设置边界类型
            ssp1d.set_Bottom_Line(5000, 1600.0, 0.2, 0.0, 0.0, 1.5);
            ssp1d.set_Top_type(BC_Mode::MODE_V_Vacuum);
            // 转换ssp
            ssp::SSPStructure ssp = convert_to_SSPStructure(ssp1d);
            params.SSP.resize(1); // 只能指针
            params.SSP[0] = ssp;
        }
    };
}
