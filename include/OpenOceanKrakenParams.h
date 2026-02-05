#ifndef PARAMS_H
#define PARAMS_H
#include <vector>
#include <string>
#include <complex>
#include <fstream>
#include <Eigen/Dense>
#include <iostream>
#include <memory>
#include <float.h>
#include "nlohmann/json.hpp"
#include "ThreadPool.h"
#include <algorithm>
#include "json_eigen.hpp"
#define _Order_json 1 // 1表示使用ordered_json，0表示使用json
#ifdef _Order_json
#define OpenOcean_json nlohmann::ordered_json
#else
#define OpenOcean_json nlohmann::json
#endif



namespace OpenOceanKraken
{
    // @brief 定义了一些常量
    const double eps = 0.00737;
    constexpr double pi = 3.14159265358979323846;
    constexpr bool ThreeD = false;
    constexpr double RadDeg = 180 / pi;
    constexpr double DegRad = pi / 180;
    constexpr double c0 = 1500;
    constexpr double HUGE1 = 1.0e8;
    constexpr int MaxSSP = 1001;
    constexpr int MaxBisections = 50;
    constexpr std::complex<float> I1(0, 1);
    constexpr std::complex<double> I1D(0.0, 1.0);
    // 定义最大迭代次数
    constexpr int MAXIT = 1;

    // 外推系数
    constexpr int NSet = 5;
    constexpr int BCIiPowerR = 50;
    constexpr int BCIiPowerF = -50;
    constexpr double BCIRoof = 1.0e+50;
    constexpr double BCIFloor = 1.0e-50;

    template <typename T>
    constexpr T SQ(const T &a) { return a * a; }

    // #ifndef USE_FLOAT
    // #define USE_FLOAT 1
    //     using REAL = double;
    //     using COMPLEX = std::complex<REAL>;
    // #else
    // #define USE_FLOAT 1
    //     using REAL = float;
    //     using COMPLEX = std::complex<REAL>;
    // #endif

#define _Order_json 1 // 1表示使用ordered_json，0表示使用json
#ifdef _Order_json
#define kkc_json nlohmann::ordered_json
#else
#define kkc_json nlohmann::json
#endif

    // ssp类型
    enum class SSP_Mode
    {
        MODE_N_n2Linear, // N2 线性插值；
        MODE_C_cLinear,  // C-线性插值；                 ——默认
        MODE_P_cPCHIP,   // 分段三次埃尔米特插值多项式
        MODE_S_cCubic,   // 三次样条插值；
        MODE_A_Analytic,
        MODE_Q_Quad, // 声速场二次逼近；要输入ssp矩阵
    };

    // ssp类型
    enum class Media_Mode
    {
        MODE_A_Acoustic, // 声学层（没有横波）
        MODE_E_Elastic,  //
    };

    // 衰减选项
    // 1. 衰减单位类型（核心物理单位）
    enum class AttenuationUnit
    {
        MODE_F_dB_per_m_kHz,   // 衰减单位采用(dB/m)kHz；
        MODE_L_params_lose,    // 衰减单位采用参数损失；
        MODE_M_dB_per_m,       // 衰减单位采用 dB/m；
        MODE_m_dB_per_m,       // 衰减单位采用 dB/m；
        MODE_N_Nepers_per_m,   // 衰减单位采用 Nepers/m；
        MODE_Q_Quality_Factor, // 衰减单位采用 Q 因子；
        MODE_W_db_per_lambda,  // 衰减单位采用 dB/λ(波长) ——默认
    };
    // 2. 海洋吸收模型（可选附加）
    enum class OceanAbsorptionModel
    {
        None,     // 普通
        Thorpe,   // T
        FrancGarr // F
    };

    struct Atten_Mode
    {
        AttenuationUnit attnUnit = AttenuationUnit::MODE_W_db_per_lambda; // 衰减单位类型
        OceanAbsorptionModel absModel = OceanAbsorptionModel::None;       // 海洋吸收模型
    };

    // 边界条件类型
    enum class BC_Mode
    {
        MODE_R_Rigid,       // 刚性边界条件
        MODE_V_Vacuum,      // 真空                       ——默认
        MODE_F_File,        // 从文件读取边界条件
        MODE_A_Half_space,  // 半空间边界条件
        MODE_G_Grain,       // 粒子边界条件
        MODE_P_Precomputed, // 预计算反射系数
    };

    // 声源类型
    enum class Source_Mode
    {
        MODE_R_Point, // 点源 ———默认
        MODE_X_Line,  // 线源
    };

    // 接收网格类型
    enum class Grid_Mode
    {
        MODE_R_Rectangular, // 矩形网格    ——默认
        MODE_I_Irregular,   // 不规则网格
    };

    // Kraken特有的运行模式
    enum class Run_Mode
    {
        MODE_M_Modes, // 只计算简正波模式
        MODE_F_Field, // 计算声场
        MODE_B_Both   // 同时计算模式和声场
    };

    // 相干和非相干
    enum class CoherenceType
    {
        Coherent,  // 相干
        Incoherent // 非相干
    };

    // ModeType
    enum class ModeType
    {
        Adiabatic, // 绝热模式
        Couple     // 耦合模式
    };

    // @brief 半空间属性结构体
    struct HSInfo
    {
        // @brief 纵波（压缩波P-wave）声速
        double alphaR;
        // @brief 纵波吸收系数
        double alphaI;
        // @brief 横波（剪切波S-wave）声速
        double betaR;
        // @brief 横波吸收系数
        double betaI;
        double beta, ft;
        // @brief P-wave速度
        std::complex<double> cp;
        // @brief S-wave速度
        std::complex<double> cs;
        // @brief 密度
        double rho;
        // @brief 深度
        double Depth;
        // @brief 边界条件类型
        BC_Mode BC;

        
    };
    // ssp结构体
    namespace ssp
    {
        // @brief 声速剖面结构体
        struct SSPStructure
        {
            // @brief 该剖面层数
            int NMedia;
            // @brief 第一个声学层索引
            int FirstAcoustic;
            // @brief 最后一个声学层索引
            int LastAcoustic;
            // 多层介质层线性存储，用NPts来指代总数据点数，NMesh指代抽样点数，NMedia来指代层数
            // 介质层定义
            Eigen::VectorXd beta;  // size = NMedia
            Eigen::VectorXd ft;    // size = NMedia
            Eigen::VectorXd sigma; // size = NMedia
            std::vector<Media_Mode> Material;
            SSP_Mode SSPType;

            // 层厚度
            Eigen::VectorXd depth;
            // @brief 声速剖面点数
            Eigen::VectorXi NPts; // size = NMedia
            // @brief多层介质存储偏移量
            Eigen::VectorXi offset; // size = NMedia
            // @brief 插值偏移量
            Eigen::VectorXi interp_offset; // size = NMedia

            // @brief 声速剖面细分点数
            Eigen::VectorXi NMesh; // size = NMedia

            // @brief 深度向量
            Eigen::VectorXd z;
            Eigen::VectorXd alphaR; // 声速，纵波速度
            Eigen::VectorXd alphaI; // 横波速度
            Eigen::VectorXd betaR;  // 纵波衰减
            Eigen::VectorXd betaI;  // 横波衰减
            Eigen::VectorXd rho;    // 密度

            // @brief 声速向量
            Eigen::VectorXcd cp; // 纵波声速，根据alphaR和alphaI算出的复声速
            Eigen::VectorXcd cs; // 横波声速，根据betaR和betaI算出的复声速

            // @brief 声速三次样条系数矩阵
            // 默认行大小
            int row_size = 4;
            // @brief 声速三次样条系数矩阵
            Eigen::MatrixXcd cspline;
            // @brief PCHIP 系数
            Eigen::MatrixXcd cCoef;
            // @brief PCHIP 系数
            Eigen::MatrixXcd csWork;
            // @brief P波PCHIP系数 (4 x MaxSSP)
            Eigen::MatrixXcd cpCoef;
            // @brief S波PCHIP系数 (4 x MaxSSP)
            Eigen::MatrixXcd csCoef;
            // @brief 密度PCHIP系数 (4 x MaxSSP)
            Eigen::MatrixXcd rhoCoef;
            // @brief P波三次样条系数
            Eigen::MatrixXcd cpSpline;
            // @brief S波三次样条系数
            Eigen::MatrixXcd csSpline;
            // @brief 密度三次样条系数
            Eigen::MatrixXcd rhoSpline;

            HSInfo HSTop; // 顶层
            HSInfo HSBot; // 底层

            // 返回第 iMedia 层在全局向量中的起始索引（包含）
            inline int get_media_start(int iMedia) const
            {
                assert(iMedia >= 0 && iMedia < NMedia);
                return offset[iMedia];
            }

            // 返回第 iMedia 层的结束索引（不包含，C++ 半开区间惯例）
            inline int get_media_end(int iMedia) const
            {
                assert(iMedia >= 0 && iMedia < NMedia);
                return offset[iMedia] + NPts[iMedia];
            }

            // 返回该层点数（冗余但方便）
            inline int get_media_size(int iMedia) const
            {
                assert(iMedia >= 0 && iMedia < NMedia);
                return NPts[iMedia];
            }
            // 返回该层网格点数
            inline int get_media_Nmesh(int iMedia) const
            {
                assert(iMedia >= 0 && iMedia < NMedia);
                return NMesh[iMedia];
            }
            inline int get_global_interp_offset(int iMedium) const
            {
                return interp_offset[iMedium];
            }
        };

        struct SSPLayer
        {
            // --- 成员变量 ---
            int npts = 0;
            int nmesh = 0;
            double beta = 0.0;
            double ft = 0.0;
            double sigma = 0.0;
            Media_Mode Material = Media_Mode::MODE_A_Acoustic;

            Eigen::VectorXd z;
            Eigen::VectorXd rho;
            Eigen::VectorXd alphaR;
            Eigen::VectorXd alphaI;
            Eigen::VectorXd betaR;
            Eigen::VectorXd betaI;

            // --- 构造函数与赋值 ---
            SSPLayer() = default;

            SSPLayer(
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
                Media_Mode media);

            SSPLayer(SSPLayer &&other) noexcept;
            SSPLayer(const SSPLayer &other) = default;

            SSPLayer &operator=(const SSPLayer &other) = default;
            SSPLayer &operator=(SSPLayer &&other) noexcept = default;

            // --- 成员函数声明 ---
            inline bool empty() const;
            inline bool is_valid() const;
            // 
        };

        // 一个距离无关区域包含SSP和底部与顶部的半空间属性
        //  ----------------------------
        //  FlattenedData 必须完整定义在头文件中！
        //  ----------------------------
        struct FlattenedData
        {
            Eigen::VectorXi NPts;
            Eigen::VectorXi NMesh;
            Eigen::VectorXd beta, ft, sigma;
            Eigen::VectorXd z, rho, alphaR, alphaI, betaR, betaI;
            // FlattenedData(FlattenedData &&) = default;
            // FlattenedData &operator=(FlattenedData &&) = default;
        };

        struct Range_Independent_Area
        {
            SSP_Mode SSPType = SSP_Mode::MODE_C_cLinear;
            std::vector<SSPLayer> layers;
            HSInfo HSTop;
            HSInfo HSBot;
            double Range = 0; //距离无关区域在声场中的位置
            // 构造函数（如有）
            Range_Independent_Area() = default;

            // 修改接口（简单，可 inline）
            inline void addLayer(const SSPLayer &layer)
            {
                layers.push_back(layer);
            }

            inline void insertLayer(size_t i, const SSPLayer &layer)
            {
                assert(i <= layers.size());
                layers.insert(layers.begin() + i, layer);
            }

            inline void removeLayer(size_t i)
            {
                assert(i < layers.size());
                layers.erase(layers.begin() + i);
            }

            inline void clearLayer()
            {
                layers.clear();
            }

            // 返回层数（冗余但方便）
            inline int NMedia() const { return static_cast<int>(layers.size()); }

            // 半空间设置（声明，实现在 .cpp）
            void set_Bottom_Line(double Depth, double alphaR, double alphaI, double betaR, double betaI, double rho); // 仅在Half-space模式下有效
            void set_Top_Line(double Depth, double alphaR, double alphaI, double betaR, double betaI, double rho);    // 仅在Half-space模式下有效
            void set_Bottom_type(BC_Mode type);
            void set_Top_type(BC_Mode type);

            // flatten() 声明
            FlattenedData flatten() const;
        };

    }

    struct Position
    {

        // @brief 声源z坐标点的数量
        int NSz;
        // @brief 接收器z坐标点的数量
        int NRz;
        // @brief 接收器r坐标点的数量
        int NRr;
        // @brief 每个距离网格的深度结构个数
        int NRz_per_range;

        // @brief 距离间隔
        double Delta_r;

        // @brief 声源点z坐标
        Eigen::VectorXd Sz;

        // @brief 接收器r坐标
        Eigen::VectorXd Rr;
        // @brief 接收器z坐标
        Eigen::VectorXd Rz;
        // @brief 阵列水平倾斜距离
        Eigen::VectorXd Ro;
        // @brief 接收水平方向角
        Eigen::VectorXd theta;
        bool is_Linspace_Rr = false; // 水平是否等间距
        bool is_Linspace_Rz = false; // 垂直是否等间距
        bool is_Linspace_Sz = false; // 声源是否等间距
        Grid_Mode GridType;
    };
    // @brief 反射系数结构体
    struct ReflectionCoef
    {
        double theta; // 掠射角
        double R;     // 反射系数幅度
        double phi;   // 反射系数相位
    };

    // @brief 声源指向性图结构体
    struct SrcBmPat
    {
        // @brief 声源指向性图点数
        int NSBPPts;
        // @brief 声源指向性图标志
        std::string SBPFlag;
        // @brief 声源指向性图向量
        Eigen::VectorXd theta;
        Eigen::VectorXd pat;
        bool isSet = false; // 是否设置了指向性图
    };

    struct FreqInfo
    {
        // @brief 频率
        double freq;
        // @brief 频率个数
        int Nfreq;
        // @brief 频率向量
        Eigen::VectorXd freqvec;
    };

    struct BdryType
    {
        // @brief 将BdryPt实例化为Top结构体
        HSInfo Top;
        // @brief 将BdryPt实例化为Bot结构体
        HSInfo Bot;
    };
    struct ReflectionCoefInfo
    {
        // @brief 顶部反射系数参数
        Eigen::Matrix<ReflectionCoef, 1, Eigen::Dynamic> RTop;
        // @brief 底部反射系数参数
        Eigen::Matrix<ReflectionCoef, 1, Eigen::Dynamic> RBot;
        bool isDeg = false; // @brief 是否是角度制
    };

    // 有限差分网格参数
    struct MeshParams
    {
        int NSets = 5;                // 网格集数量
        int NV[5] = {1, 2, 4, 8, 16}; // Richardson外推系数数组
    };

    // 本征值相关参数
    struct EigenParams
    {
        int M;                    // 模式数量
        int firstM;               // iset=0时的模式个数，决定了矩阵维度大小 2 D f / cmin * 1.1
        Eigen::VectorXd EVMat;    // 本征值矩阵(一维向量化)
        Eigen::VectorXd Extrap;   // 外推矩阵
        Eigen::VectorXcd k;       // 波数向量
        Eigen::VectorXd VG;       // 群速度向量
        Eigen::MatrixXcd PsiR;    // 接收器深度本征函数值
        Eigen::MatrixXcd PsiS;    // 声源深度本征函数值
        Eigen::MatrixXcd dPsidzR; // 接收器深度本征函数值对深度微分
        Eigen::MatrixXcd dPsidzS; // 声源深度本征函数值对深度微分
        // VectorXi modes;    // 模式索引

        inline void resize(int firstM, int NSz, int NRz, int NMeshMax, int NSets)
        {
            EVMat.resize(firstM * NSets);
            Extrap.resize(firstM * NSets);
            k.resize(firstM);
            VG.resize(firstM);
            PsiR.resize(firstM, NRz);
            PsiS.resize(firstM, NSz);
            dPsidzR.resize(firstM, NRz);
            dPsidzS.resize(firstM, NSz);
        }

        inline void setZero()
        {
            EVMat.setZero();
            Extrap.setZero();
            k.setZero();
            VG.setZero();
            PsiR.setZero();
            PsiS.setZero();
            dPsidzR.setZero();
            dPsidzS.setZero();
        }
    };

    struct OOK_parameters
    {
        // @brief 标题
        std::string Title;

        // @brief 任务数（声源*发射声线个数）
        int totalTasks;

        // @brief 距离剖面向量
        Eigen::VectorXd RProf;
        int NMeshMax;  // 最大网格数
        int NMediaMax; // 最大媒质数

        // @brief 频率信息
        FreqInfo freqinfo;

        // // @brief 媒质数
        // int NMedia;

        // @brief 吸收单位
        Atten_Mode AttenUnit;

        // @brief 声源位置
        Position Pos;

        // @brief 声速剖面参数
        std::vector<ssp::SSPStructure> SSP;
        std::vector<ssp::Range_Independent_Area> sspInput; // SSP输入 方便to_json
        // @brief 边界参数
        BdryType Bdry;

        ReflectionCoefInfo ReflectionCoef; // @brief 反射系数参数

        // @brief 实例化SrcBmPat为SBP
        SrcBmPat SBP;

        bool is_Velocity = false; // @brief 是否计算振速
        // kkc_Log *log;             // 日志

        // 核心参数
        double cLow;            // 最小相速度
        double cHigh;           // 最大相速度
        double Rmax;            // 最大距离
        Source_Mode SourceType; // 声源模式

        // 网格参数
        MeshParams mesh; // 网格参数

        // 运行模式
        Run_Mode runMode; // Kraken运行模式

        // 相干和非相干
        CoherenceType coherenceType; // 相干和非相干

        // 模式类型
        ModeType modeType; // 模式类型

    };

    // @brief 三对角矩阵结构体，kraken计算重要的中间变量
    struct TridMtx
    {
        // @brief 深度向量
        Eigen::VectorXd z;
        // @brief 声速剖面插值点声速
        Eigen::VectorXcd cp_int;
        Eigen::VectorXcd cs_int;
        Eigen::VectorXd rho_int;

        // 计算矩阵
        Eigen::VectorXd B1;
        Eigen::VectorXd B1C;
        Eigen::VectorXd B2;
        Eigen::VectorXd B3;
        Eigen::VectorXd B4;
        Eigen::VectorXd rho;

        // 本征函数
        Eigen::VectorXd psi;    // 原始mesh的本征函数值
        Eigen::VectorXd dpsidz; // 原始mesh的本征函数值的z方向导数

        // @brief 网格参数
        Eigen::VectorXi N;   // 各层的网格点数
        Eigen::VectorXd h;   // 各层的网格步长
        Eigen::VectorXd hV;  // 网格步长向量
        Eigen::VectorXi Loc; // 各层的网格点索引

        // @brief 相速度范围
        double cLow;
        double cHigh;

        inline void resize(int Maxsize, int NMediaMax, int NSets)
        {
            z.resize(Maxsize);
            cp_int.resize(Maxsize);
            cs_int.resize(Maxsize);
            rho_int.resize(Maxsize);
            B1.resize(Maxsize);
            B1C.resize(Maxsize);
            B2.resize(Maxsize);
            B3.resize(Maxsize);
            B4.resize(Maxsize);
            rho.resize(Maxsize);
            psi.resize(Maxsize);
            dpsidz.resize(Maxsize);
            N.resize(NMediaMax);
            h.resize(NMediaMax);
            hV.resize(NSets);
            Loc.resize(NMediaMax);
        }
    };

    // 输出结构
    struct OOK_output
    {
        // 本征值参数
        EigenParams *eigen;                // 本征值参数
        std::complex<float> *u_AllSources; // 声压，一维化存储，内存优化较好
        std::complex<float> *v_AllSources; // 垂直振速，一维化存储，内存优化较好
        std::complex<float> *h_AllSources; // 水平振速，一维化存储，内存优化较好
    };
}

#endif // PARAMS_H
