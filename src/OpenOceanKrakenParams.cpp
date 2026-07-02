// 结构体函数实现
#include "OpenOceanKrakenParams.h"

namespace OpenOceanKraken
{

    namespace ssp
    {

        // @brief 声速层结构体
        SSPLayer::SSPLayer(
            int npoints,
            int nmesh,
            double beta,
            double ft,
            double sigma,
            const Eigen::Ref<const Eigen::VectorXd> &_z,
            const Eigen::Ref<const Eigen::VectorXd> &_aR,
            const Eigen::Ref<const Eigen::VectorXd> &_aI,
            const Eigen::Ref<const Eigen::VectorXd> &_bR,
            const Eigen::Ref<const Eigen::VectorXd> &_bI,
            const Eigen::Ref<const Eigen::VectorXd> &_rho,
            Media_Mode media)
            : npts(npoints), nmesh(nmesh), beta(beta), ft(ft), sigma(sigma), z(_z), rho(_rho), alphaR(_aR), alphaI(_aI), betaR(_bR), betaI(_bI), Material(media)
        {
            assert(npts >= 0 && "npts must be non-negative");
            assert(nmesh >= 0 && "nmesh must be non-negative");
            assert(z.size() == npts);
            assert(rho.size() == npts);
            assert(alphaR.size() == npts);
            assert(alphaI.size() == npts);
            assert(betaR.size() == npts);
            assert(betaI.size() == npts);
        }

        // 移动构造函数
        SSPLayer::SSPLayer(SSPLayer &&other) noexcept
            : npts(other.npts), nmesh(other.nmesh), beta(other.beta), ft(other.ft), sigma(other.sigma), z(std::move(other.z)), rho(std::move(other.rho)), alphaR(std::move(other.alphaR)), alphaI(std::move(other.alphaI)), betaR(std::move(other.betaR)), betaI(std::move(other.betaI)), Material(other.Material)
        {
            other.npts = 0;
            other.nmesh = 0;
        }

        // 辅助函数
        inline bool SSPLayer::empty() const
        {
            return npts == 0;
        }

        inline bool SSPLayer::is_valid() const
        {
            return npts >= 0 &&
                   nmesh >= 0 &&
                   z.size() == npts &&
                   rho.size() == npts &&
                   alphaR.size() == npts &&
                   alphaI.size() == npts &&
                   betaR.size() == npts &&
                   betaI.size() == npts;
        }

        void Range_Independent_Area::set_Bottom_Line(double Depth, double alphaR, double alphaI, double betaR, double betaI, double rho)
        {
            this->HSBot.Depth = Depth;
            this->HSBot.alphaR = alphaR;
            this->HSBot.alphaI = alphaI;
            this->HSBot.betaR = betaR;
            this->HSBot.betaI = betaI;
            this->HSBot.rho = rho;
            //    this->HSBot.sigma = 0.0;
            this->HSBot.ft = 0.0;
            this->HSBot.BC = BC_Mode::MODE_A_Half_space;
        }

        void Range_Independent_Area::set_Top_Line(double Depth, double alphaR, double alphaI, double betaR, double betaI, double rho)
        {
            this->HSTop.Depth = Depth;
            this->HSTop.alphaR = alphaR;
            this->HSTop.alphaI = alphaI;
            this->HSTop.betaR = betaR;
            this->HSTop.betaI = betaI;
            this->HSTop.rho = rho;
            this->HSTop.ft = 0.0;
            this->HSTop.BC = BC_Mode::MODE_A_Half_space;
        }
        void Range_Independent_Area::set_Bottom_type(BC_Mode type)
        {
            this->HSBot.BC = type;
        }
        void Range_Independent_Area::set_Top_type(BC_Mode type)
        {
            this->HSTop.BC = type;
        }


        FlattenedData Range_Independent_Area::flatten() const
        {
            if (layers.empty())
            {
                return {};
            }

            int nmedia = static_cast<int>(layers.size());
            FlattenedData flat;

            // 1. 层元数据
            flat.NPts.resize(nmedia);
            flat.NMesh.resize(nmedia);
            flat.beta.resize(nmedia);
            flat.ft.resize(nmedia);
            flat.sigma.resize(nmedia);

            int total_pts = 0;
            for (int i = 0; i < nmedia; ++i)
            {
                flat.NPts[i] = layers[i].npts;
                flat.NMesh[i] = layers[i].nmesh;
                flat.beta[i] = layers[i].beta;
                flat.ft[i] = layers[i].ft;
                flat.sigma[i] = layers[i].sigma;
                total_pts += layers[i].npts;
            }

            // 2. 全局物理量
            flat.z.resize(total_pts);
            flat.rho.resize(total_pts);
            flat.alphaR.resize(total_pts);
            flat.alphaI.resize(total_pts);
            flat.betaR.resize(total_pts);
            flat.betaI.resize(total_pts);

            int offset = 0;
            for (const auto &layer : layers)
            {
                flat.z.segment(offset, layer.npts) = layer.z;
                flat.rho.segment(offset, layer.npts) = layer.rho;
                flat.alphaR.segment(offset, layer.npts) = layer.alphaR;
                flat.alphaI.segment(offset, layer.npts) = layer.alphaI;
                flat.betaR.segment(offset, layer.npts) = layer.betaR;
                flat.betaI.segment(offset, layer.npts) = layer.betaI;
                offset += layer.npts;
            }

            return flat;
        }
    }

}
