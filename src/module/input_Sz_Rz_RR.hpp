#include "paramsBase.h"

class input_Sz_Rz_RR : public paramsBase
{
public:
    input_Sz_Rz_RR() {}
    virtual ~input_Sz_Rz_RR() {}

    virtual void Init(parameters &params) const override
    {
    }

    virtual void Default(parameters &params) const override
    {
        setMunk(params);
        // std::cout<<"input_Sz_Rz_RR::Default() is called"<<std::endl;
    }

    virtual void Preprocess(parameters &params) const override
    {
        auto &pos = params.Pos;
        pos->Delta_r = 0.0;
        if (pos->NRr > 2)
        {
            pos->Delta_r = (pos->Rr[pos->NRr - 1] - pos->Rr[pos->NRr - 2]);
        }

        pos->NRz_per_range = (pos->GridType == Grid_Mode::MODE_I_Irregular) ? 1 : pos->NRz;

        // std::cout<<"input_Sz_Rz_RR::Preprocess() is called"<<std::endl;
    }

    // 手动设置全部水平接收
    void set_Rr(parameters &params, const VectorXd &Rr)
    {
        auto &pos = params.Pos;
        pos->NRr = Rr.size();
        pos->Rr.resize(pos->NRr);
        pos->Rr = Rr;
        pos->is_Linspace_Rr = false;
    }

    // 插值生成水平接收
    void set_Rr(parameters &params, const double &start, const double &end, const int &NRr)
    {
        auto &pos = params.Pos;
        pos->NRr = NRr;
        pos->Rr.resize(pos->NRr);
        linspace(start, end, pos->NRr, pos->Rr); // 线性插值
        pos->is_Linspace_Rr = true;
    }

    // 手动设置全部垂直接收
    void set_Rz(parameters &params, const VectorXd &Rz)
    {
        auto &pos = params.Pos;
        pos->NRz = Rz.size();
        pos->Rz.resize(pos->NRz);
        pos->Rz = Rz;
        pos->is_Linspace_Rz = false;
    }

    // 插值生成垂直接收
    void set_Rz(parameters &params, const double &start, const double &end, const int &NRz)
    {
        auto &pos = params.Pos;
        pos->NRz = NRz;
        pos->Rz.resize(pos->NRz);
        linspace(start, end, pos->NRz, pos->Rz); // 线性插值
        pos->is_Linspace_Rz = true;
    }

    // 手动设置声源
    void set_Sz(parameters &params, const VectorXd &Sz)
    {
        auto &pos = params.Pos;
        pos->NSz = Sz.size();
        pos->Sz.resize(pos->NSz);
        pos->Sz = Sz;
        pos->is_Linspace_Sz = false;
    }
    // 插值生成声源
    void set_Sz(parameters &params, const double &start, const double &end, const int &NSz)
    {
        auto &pos = params.Pos;
        pos->NSz = NSz;
        pos->Sz.resize(pos->NSz);
        linspace(start, end, pos->NSz, pos->Sz); // 线性插值
        pos->is_Linspace_Sz = true;
    }

    // 设置网格类型
    void set_GridType(parameters &params, const Grid_Mode &GridType)
    {
        auto &pos = params.Pos;
        pos->GridType = GridType;
    }

    VectorXd get_Rr(const parameters &params)
    {
        auto &pos = params.Pos;
        return pos->Rr;
    }

    VectorXd get_Rz(const parameters &params)
    {
        auto &pos = params.Pos;
        return pos->Rz;
    }

    VectorXd get_Sz(const parameters &params)
    {
        auto &pos = params.Pos;
        return pos->Sz;
    }

    // 获取网格类型
    Grid_Mode get_GridType(const parameters &params)
    {
        auto &pos = params.Pos;
        return pos->GridType;
    }

    // 打印
    virtual void Echo(parameters &params) const override
    {
    }

private:
    // 设置默认值的私有方法
    void setDefaultValues(parameters &params) const
    {
        auto &pos = params.Pos;
        pos->NSz = 1;
        pos->NSx = 1;
        pos->NSy = 1;
        pos->NRz = 201;
        pos->NRr = 300;

        pos->Sz.resize(pos->NSz);
        pos->Sz(0) = 25;

        // 默认水平接收
        pos->Rr.resize(pos->NRr);
        for (int32_t i = 0; i < pos->NRr; ++i)
        {
            pos->Rr[i] = 100.0 * (double)(i + 1);
        }

        pos->Rz.resize(pos->NRz);
        pos->Ro.resize(pos->NRz);
        for (int32_t i = 0; i < pos->NRz; ++i)
        {
            pos->Rz[i] = i;
            pos->Ro[i] = 0;
        }
        pos->GridType = Grid_Mode::MODE_R_Rectangular;

        pos->is_Linspace_Rr = false;
        pos->is_Linspace_Rz = false;
        pos->is_Linspace_Sz = false;
    }

    void setMunk(parameters &params) const
    {
        // 声源接收设置
        auto &pos = params.Pos;
        pos->NSz = 1; // 对应env中的NSD=1
        pos->Sz.resize(pos->NSz);
        pos->Sz(0) = 1000.0; // 对应env中的SD(1:NSD)=1000.0 (m)

        pos->NRr = 501; // 对应env中的NR=1001
        pos->Rr.resize(pos->NRr);
        for (size_t i = 0; i < pos->NRr; ++i)
        {
            // 对应env中的R范围0.0-100.0 km，转换为m并线性分布
            pos->Rr(i) = 100.0 * i; // 步长100m (100000m / 1000步)
        }

        pos->NRz = 1001; // 对应env中的NRD=501
        pos->NRz_per_range = pos->NRz;
        pos->Rz.resize(pos->NRz);
        pos->Ro.resize(pos->NRz);
        for (size_t i = 0; i < pos->NRz; ++i)
        {
            // 对应env中的RD范围0.0-5000.0 m，线性分布
            pos->Rz(i) = 5.0 * (i); // 步长10m (5000m / 500步)
            pos->Ro(i) = 0;
        }
        pos->GridType = Grid_Mode::MODE_R_Rectangular;

        pos->is_Linspace_Rr = false;
        pos->is_Linspace_Rz = false;
        pos->is_Linspace_Sz = false;
    }
};
