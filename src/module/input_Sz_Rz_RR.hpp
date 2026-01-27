#include "paramsBase.h"
namespace OpenOceanKraken
{
    class input_Sz_Rz_RR : public paramsBase
    {
    public:
        input_Sz_Rz_RR() {}
        virtual ~input_Sz_Rz_RR() {}

        virtual void Init(OOK_parameters &params) const override
        {
        }

        virtual void Default(OOK_parameters &params) const override
        {
            setMunk(params);
            // std::cout<<"input_Sz_Rz_RR::Default() is called"<<std::endl;
        }

        virtual void Preprocess(OOK_parameters &params) const override
        {
            auto &pos = params.Pos;
            pos.Delta_r = 0.0;
            if (pos.NRr > 2)
            {
                pos.Delta_r = (pos.Rr[pos.NRr - 1] - pos.Rr[pos.NRr - 2]);
            }

            pos.NRz_per_range = (pos.GridType == Grid_Mode::MODE_I_Irregular) ? 1 : pos.NRz;

            // std::cout<<"input_Sz_Rz_RR::Preprocess() is called"<<std::endl;
        }

        // 手动设置全部水平接收
        void set_Rr(OOK_parameters &params, const Eigen::VectorXd &Rr)
        {
            auto &pos = params.Pos;
            pos.NRr = Rr.size();
            pos.Rr.resize(pos.NRr);
            pos.Rr = Rr;
            pos.is_Linspace_Rr = false;
        }

        // 插值生成水平接收
        void set_Rr(OOK_parameters &params, const double &start,const double &end, const int &NRr)
        {
            auto &pos = params.Pos;
            pos.NRr = NRr;
            pos.Rr.resize(pos.NRr);
            Util::linspace(start, end, pos.NRr, pos.Rr); // 线性插值
            pos.is_Linspace_Rr = true;
        }

        // 手动设置全部垂直接收
        void set_Rz(OOK_parameters &params, const Eigen::VectorXd &Rz)
        {
            auto &pos = params.Pos;
            pos.NRz = Rz.size();
            pos.Rz.resize(pos.NRz);
            pos.Rz = Rz;
            pos.is_Linspace_Rz = false;
        }

        // 插值生成垂直接收
        void set_Rz(OOK_parameters &params, const double &start, const double &end, const int &NRz)
        {
            auto &pos = params.Pos;
            pos.NRz = NRz;
            pos.Rz.resize(pos.NRz);
            Util::linspace(start, end, pos.NRz, pos.Rz); // 线性插值
            pos.is_Linspace_Rz = true;
        }

        // 手动设置声源
        void set_Sz(OOK_parameters &params, const Eigen::VectorXd &Sz)
        {
            auto &pos = params.Pos;
            pos.NSz = Sz.size();
            pos.Sz.resize(pos.NSz);
            pos.Sz = Sz;
            pos.is_Linspace_Sz = false;
        }
        // 插值生成声源
        void set_Sz(OOK_parameters &params, const double &start, const double &end, const int &NSz)
        {
            auto &pos = params.Pos;
            pos.NSz = NSz;
            pos.Sz.resize(pos.NSz);
            Util::linspace(start, end, pos.NSz, pos.Sz); // 线性插值
            pos.is_Linspace_Sz = true;
        }

        // 设置网格类型
        void set_GridType(OOK_parameters &params, const Grid_Mode &GridType)
        {
            auto &pos = params.Pos;
            pos.GridType = GridType;
        }

        // 打印
        virtual void Echo(OOK_parameters &params) const override
        {
        }

    private:
        // 设置默认值的私有方法

        void setMunk(OOK_parameters &params) const
        {
            auto &pos = params.Pos;

            // === 声源设置 ===
            constexpr int NSz_default = 1;
            pos.NSz = NSz_default;
            pos.Sz.resize(NSz_default);
            pos.Sz(0) = 1000.0;         // 深度 1000 米
            pos.is_Linspace_Sz = false; // 单点，非 linspace

            // === 水平接收距离 (Range) ===
            constexpr int NRr_default = 601;
            constexpr double R_start = 0.0;   // 起始距离 (m)
            constexpr double R_end = 100.0e3; // 100 km -> 米

            pos.NRr = NRr_default;
            pos.Rr.resize(NRr_default);
            Util::linspace(R_start, R_end, NRr_default, pos.Rr);
            pos.is_Linspace_Rr = true; // ✅ 修正：确实是线性分布！

            // === 垂直接收深度 (Depth) ===
            constexpr int NRz_default = 1001;
            constexpr double Z_max = 5000.0; // 最大深度 5000 米

            pos.NRz = NRz_default;
            pos.NRz_per_range = NRz_default;
            pos.Rz.resize(NRz_default);
            pos.Ro.resize(NRz_default);

            // 使用 linspace 保持风格一致（更安全、简洁）
            Util::linspace(0.0, Z_max, NRz_default, pos.Rz);

            // Ro 全零（方位角？）
            pos.Ro.setZero(); // 比循环赋值更高效

            pos.is_Linspace_Rz = true; // ✅ 修正：确实是线性分布！

            // === 网格类型 ===
            pos.GridType = Grid_Mode::MODE_R_Rectangular;
        }
    };
}