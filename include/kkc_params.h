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

#define SQ(a) ((a) * (a)) // Square
using namespace std;
using namespace Eigen;
// using namespace nlohmann;

#define _Order_json 1 // 1表示使用ordered_json，0表示使用json
#ifdef _Order_json
#define kkc_json nlohmann::ordered_json
#else
#define kkc_json nlohmann::json
#endif

class kkc_Log; // 前向声明

// @brief 定义了一些常量
const double eps = 0.00737;
constexpr double pi = 3.14159265358979323846;
constexpr bool ThreeD = false;
constexpr double RadDeg = 180 / pi;
constexpr double DegRad = pi / 180;
constexpr double c0 = 1500;
constexpr double HUGE1 = 1.0e8;
constexpr int MaxSSP = 1001;
const int MaxBisections = 50;
constexpr std::complex<float> I1(0, 1);
constexpr std::complex<double> I1D(0.0, 1.0);
// 定义最大迭代次数
const int MAXIT = 1;

// 外推系数
const int NSet = 5;
const int BCIiPowerR = 50;
const int BCIiPowerF = -50;
const double BCIRoof = 1.0e+50;
const double BCIFloor = 1.0e-50;

// 选项的枚举

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
enum class Atten_Mode
{
    MODE_F_dB_per_m_kHz,   // 衰减单位采用(dB/m)kHz；
    MODE_L_params_lose,    // 衰减单位采用参数损失；
    MODE_M_dB_per_m,       // 衰减单位采用 dB/m；
    MODE_m_dB_per_m,       // 衰减单位采用 dB/m；
    MODE_N_Nepers_per_m,   // 衰减单位采用 Nepers/m；
    MODE_Q_Quality_Factor, // 衰减单位采用 Q 因子；
    MODE_W_db_per_lambda,  // 衰减单位采用 dB/λ(波长) ——默认

    // 描述水体Thorpe体积衰减系数（用的比较少）
    MODE_FT_dB_per_m_kHz,   // 衰减单位采用(dB/m)kHz；
    MODE_LT_params_lose,    // 衰减单位采用参数损失；
    MODE_MT_dB_per_m,       // 衰减单位采用 dB/m；
    MODE_mT_dB_per_m,       // 衰减单位采用 dB/m；
    MODE_NT_Nepers_per_m,   // 衰减单位采用 Nepers/m；
    MODE_QT_Quality_Factor, // 衰减单位采用 Q 因子；
    MODE_WT_db_per_lambda,  // 衰减单位采用 dB/λ(波长)

    // Franc_Garr（用的比较少）
    MODE_FF_dB_per_m_kHz,   // 衰减单位采用(dB/m)kHz；
    MODE_LF_params_lose,    // 衰减单位采用参数损失；
    MODE_MF_dB_per_m,       // 衰减单位采用 dB/m；
    MODE_mF_dB_per_m,       // 衰减单位采用 dB/m；
    MODE_NF_Nepers_per_m,   // 衰减单位采用 Nepers/m；
    MODE_QF_Quality_Factor, // 衰减单位采用 Q 因子；
    MODE_WF_db_per_lambda,  // 衰减单位采用 dB/λ(波长)
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

struct rxyz_vector
{
    VectorXd r;
    VectorXd x;
    VectorXd y;
    VectorXd z;
};

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
    VectorXd beta;  // size = NMedia
    VectorXd ft;    // size = NMedia
    VectorXd sigma; // size = NMedia
    std::vector<Media_Mode> Material;
    SSP_Mode SSPType;

    // 层厚度
    VectorXd depth;
    // @brief 声速剖面点数
    VectorXi NPts; // size = NMedia
    // @brief多层介质存储偏移量
    VectorXi offset; // size = NMedia

    // @brief 声速剖面细分点数
    VectorXi NMesh; // size = NMedia

    // @brief 深度向量
    VectorXd z;
    VectorXd alphaR; // 声速，纵波速度
    VectorXd alphaI; // 横波速度
    VectorXd betaR;  // 纵波衰减
    VectorXd betaI;  // 横波衰减
    VectorXd rho;    // 密度

    // @brief 声速向量
    VectorXcd cp; // 纵波声速，根据alphaR和alphaI算出的复声速
    VectorXcd cs; // 横波声速，根据betaR和betaI算出的复声速

    // @brief 声速三次样条系数矩阵
    // 默认行大小
    int row_size = 4;
    // @brief 声速三次样条系数矩阵
    MatrixXcd cspline;
    // @brief PCHIP 系数
    MatrixXcd cCoef;
    // @brief PCHIP 系数
    MatrixXcd csWork;
    // @brief P波PCHIP系数 (4 x MaxSSP)
    MatrixXcd cpCoef;
    // @brief S波PCHIP系数 (4 x MaxSSP)
    MatrixXcd csCoef;
    // @brief 密度PCHIP系数 (4 x MaxSSP)
    MatrixXcd rhoCoef;
    // @brief P波三次样条系数
    MatrixXcd cpSpline;
    // @brief S波三次样条系数
    MatrixXcd csSpline;
    // @brief 密度三次样条系数
    MatrixXcd rhoSpline;

    // 返回第 iMedia 层在全局向量中的起始索引（包含）
    int get_media_start(int iMedia) const
    {
        assert(iMedia >= 0 && iMedia < NMedia);
        return offset[iMedia];
    }

    // 返回第 iMedia 层的结束索引（不包含，C++ 半开区间惯例）
    int get_media_end(int iMedia) const
    {
        assert(iMedia >= 0 && iMedia < NMedia);
        return offset[iMedia] + NPts[iMedia];
    }

    // 返回该层点数（冗余但方便）
    int get_media_size(int iMedia) const
    {
        assert(iMedia >= 0 && iMedia < NMedia);
        return NPts[iMedia];
    }
    // 返回该层网格点数
    int get_media_Nmesh(int iMedia) const
    {
        assert(iMedia >= 0 && iMedia < NMedia);
        return NMesh[iMedia];
    }
    int get_global_interp_offset(int iMedium) const
    {
        int offset = 0;
        for (int i = 0; i < iMedium; ++i)
        {
            offset += NMesh[i] + 1; // 每层有 Nmesh+1 个插值点
        }
        return offset;
    }
};


// 输入的SSP
struct SSPLayer {
    // --- 成员变量 ---
    int npts = 0;           // 该层原始数据点数
    int nmesh = 0;          // 该层网格点数（插值用）
    double beta = 0.0;
    double ft = 0.0;
    double sigma = 0.0;
    Media_Mode Material = Media_Mode::MODE_A_Acoustic;
    // 物理量向量（长度应等于 npts）
    Eigen::VectorXd z;
    Eigen::VectorXd rho;
    Eigen::VectorXd alphaR; // 纵波速度实部
    Eigen::VectorXd alphaI; // 纵波速度虚部（或衰减）
    Eigen::VectorXd betaR;  // 横波速度实部
    Eigen::VectorXd betaI;  // 横波速度虚部

    // --- 默认构造函数 ---
    SSPLayer() = default;

    // --- 主构造函数（推荐使用）---
    SSPLayer(
        int np,
        int nm,
        double b,
        double f,
        double s,
        const Eigen::Ref<const Eigen::VectorXd>& _z,
        const Eigen::Ref<const Eigen::VectorXd>& _rho,
        const Eigen::Ref<const Eigen::VectorXd>& _aR,
        const Eigen::Ref<const Eigen::VectorXd>& _aI,
        const Eigen::Ref<const Eigen::VectorXd>& _bR,
        const Eigen::Ref<const Eigen::VectorXd>& _bI,
        Media_Mode media
    )
        : npts(np)
        , nmesh(nm)
        , beta(b)
        , ft(f)
        , sigma(s)
        , z(_z)
        , rho(_rho)
        , alphaR(_aR)
        , alphaI(_aI)
        , betaR(_bR)
        , betaI(_bI)
        , Material(media)
    {
        // 断言：所有向量长度必须等于 npts
        assert(npts >= 0 && "npts must be non-negative");
        assert(nmesh >= 0 && "nmesh must be non-negative");
        assert(z.size() == npts);
        assert(rho.size() == npts);
        assert(alphaR.size() == npts);
        assert(alphaI.size() == npts);
        assert(betaR.size() == npts);
        assert(betaI.size() == npts);
    }

    // --- 移动构造函数（可选，但推荐）---
    SSPLayer(SSPLayer&& other) noexcept
        : npts(other.npts)
        , nmesh(other.nmesh)
        , beta(other.beta)
        , ft(other.ft)
        , sigma(other.sigma)
        , z(std::move(other.z))
        , rho(std::move(other.rho))
        , alphaR(std::move(other.alphaR))
        , alphaI(std::move(other.alphaI))
        , betaR(std::move(other.betaR))
        , betaI(std::move(other.betaI))
        , Material(other.Material)
    {
        // 可选：将 other 置为有效但空状态
        other.npts = 0;
        other.nmesh = 0;
    }

    // --- 拷贝构造函数（默认即可，但显式声明更清晰）---
    SSPLayer(const SSPLayer& other) = default;

    // --- 赋值运算符 ---
    SSPLayer& operator=(const SSPLayer& other) = default;
    SSPLayer& operator=(SSPLayer&& other) noexcept = default;

    // --- 辅助函数：检查是否为空 ---
    bool empty() const {
        return npts == 0;
    }

    // --- 辅助函数：验证内部一致性（可用于调试）---
    bool is_valid() const {
        return
            npts >= 0 &&
            nmesh >= 0 &&
            z.size() == npts &&
            rho.size() == npts &&
            alphaR.size() == npts &&
            alphaI.size() == npts &&
            betaR.size() == npts &&
            betaI.size() == npts;
    }
};


struct SSP_1D {
    SSP_Mode SSPType = SSP_Mode::MODE_C_cLinear;
    std::vector<SSPLayer> layers;  // 核心存储！
    

    void add(const SSPLayer& layer) {
        layers.push_back(layer);
    }

    void insert(size_t i, const SSPLayer& layer) {
        assert(i <= layers.size());
        layers.insert(layers.begin() + i, layer);
    }

    void remove(size_t i) {
        assert(i < layers.size());
        layers.erase(layers.begin() + i);
    }

    void clear() {
        layers.clear();
    }

    int NMedia() const { return static_cast<int>(layers.size()); }

    // ----------------------------
    // 惰性 flatten：按需生成连续向量（用于计算）
    // ----------------------------
    struct FlattenedData {
        Eigen::VectorXi NPts;
        Eigen::VectorXi NMesh;
        Eigen::VectorXd beta, ft, sigma;

        Eigen::VectorXd z, rho, alphaR, alphaI, betaR, betaI;
    };

    FlattenedData flatten() const {
        if (layers.empty()) return {};

        int nmedia = layers.size();
        FlattenedData flat;
        
        // 1. 层元数据
        flat.NPts.resize(nmedia);
        flat.NMesh.resize(nmedia);
        flat.beta.resize(nmedia);
        flat.ft.resize(nmedia);
        flat.sigma.resize(nmedia);

        int total_pts = 0;
        for (int i = 0; i < nmedia; ++i) {
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
        for (const auto& layer : layers) {
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

    // ----------------------------
    // 如果旧接口需要直接访问连续向量（兼容 legacy code）
    // 可提供只读视图（但不推荐频繁调用）
    // ----------------------------
    Eigen::Map<const Eigen::VectorXd> view_z() const {
        static thread_local Eigen::VectorXd cache; // 或由外部管理
        cache = flatten().z;
        return Eigen::Map<const Eigen::VectorXd>(cache.data(), cache.size());
    }
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
    double Mz; // 设置Grain size才需要
    // @brief 界面粗糙度
    double sigma;
};

struct BdryPt2
{
    // @brief
    HSInfo HS;
};

struct BdryType
{
    // @brief 将BdryPt实例化为Top结构体
    BdryPt2 Top;
    // @brief 将BdryPt实例化为Bot结构体
    BdryPt2 Bot;
};

// @brief 将BdryType实例化为Bdry结构体，Bdry.Top.HS.alphaR，四层结构
// BdryType Bdry;

// @brief 收发位置的结构体
struct Position
{
    // @brief 声源x坐标点的数量
    int NSx;
    // @brief 声源y坐标点的数量
    int NSy;
    // @brief 声源z坐标点的数量
    int NSz;
    // @brief 接收器z坐标点的数量
    int NRz;
    // @brief 接收器r坐标点的数量
    int NRr;
    // @brief 声源辐射水平方向角数量
    int Ntheta;
    // @brief 每个距离网格的深度结构个数
    int NRz_per_range;

    // @brief 距离间隔
    double Delta_r;
    // @brief 水平方向角度间隔
    double Delta_theta;

    // VectorXi iSz;
    // VectorXi iRz;
    // @brief 声源点x坐标
    // VectorXd Sx;
    // // @brief 声源点y坐标
    // VectorXd Sy;
    // @brief 声源点z坐标
    VectorXd Sz;

    // @brief 接收器r坐标
    VectorXd Rr;
    // @brief 接收器z坐标
    VectorXd Rz;
    // @brief 阵列水平倾斜距离
    VectorXd Ro;
    // // @brief 用于插值的权重ws
    // VectorXd ws;
    // // @brief 用于插值的权重wr
    // VectorXd wr;
    // @brief 接收水平方向角
    VectorXd theta;
    bool is_Linspace_Rr = false; // 水平是否等间距
    bool is_Linspace_Rz = false; // 垂直是否等间距
    bool is_Linspace_Sz = false; // 声源是否等间距
    Grid_Mode GridType;
};

// @brief Struct representing compressional and shear wave speeds/attenuations in user units.
// @brief 用户单位下的纵波和横波速度/衰减结构体
struct HSInfo2
{
    // @brief 纵波（压缩波P-wave）声速
    double alphaR;
    // @brief 纵波吸收系数
    double alphaI;
    // @brief 横波（剪切波S-wave）声速
    double betaR;
    // @brief 横波吸收系数
    double betaI;
    // @brief P-wave速度
    std::complex<double> cp;
    // @brief S-wave速度
    std::complex<double> cs;
    // @brief 密度
    double rho;
    // @brief 深度
    double Depth;
    // @brief 边界条件类型
    char BC;
};

// @brief 边界形状结构体
struct BdryPt
{
    // @brief 线段的坐标
    Vector2d x;
    // @brief 线段的切线
    Vector2d t;
    // @brief 线段的外法线
    Vector2d n;
    // @brief 节点处的切线（如果使用曲线坐标选项）
    Vector2d Nodet;
    // @brief 节点处的法线（如果使用曲线坐标选项）
    Vector2d Noden;
    // @brief 线段的长度
    double Len;
    // @brief 线段的曲率
    double Kappa;
    // @brief 深度的一阶导数
    double Dx;
    // @brief 深度的二阶导数
    double Dxx;
    // @brief 沿切线方向的二阶导数
    double Dss;
    // @brief 实例化HSInfo2为HS
    HSInfo2 HS;
};

// @brief 实例化BdryPt为Top和Bot可变数组
// std::vector<BdryPt> Top, Bot;

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
    VectorXd theta;
    VectorXd pat;
    bool isSet; // 是否设置了指向性图
};

struct FreqInfo
{
    // @brief 频率
    double freq;
    // @brief 频率个数
    int Nfreq;
    // @brief 频率向量
    VectorXd freqvec;
};

struct ReflectionCoefInfo
{
    // @brief 顶部反射系数参数
    Matrix<ReflectionCoef, 1, Dynamic> RTop;
    // @brief 底部反射系数参数
    Matrix<ReflectionCoef, 1, Dynamic> RBot;
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
    int M;             // 模式数量
    int firstM;        // iset=0时的模式个数，决定了矩阵维度大小 2 D f / cmin * 1.1
    VectorXd EVMat;    // 本征值矩阵(一维向量化)
    VectorXd Extrap;   // 外推矩阵
    VectorXcd k;       // 波数向量
    VectorXd VG;       // 群速度向量
    MatrixXcd PsiR;    // 接收器深度本征函数值
    MatrixXcd PsiS;    // 声源深度本征函数值
    MatrixXcd dPsidzR; // 接收器深度本征函数值对深度微分
    MatrixXcd dPsidzS; // 声源深度本征函数值对深度微分
    // VectorXi modes;    // 模式索引

    void resize(int firstM, int NSz, int NRz, int NMeshMax, int NSets)
    {
        EVMat.resize(firstM*NSets);
        Extrap.resize(firstM*NSets);
        k.resize(firstM);
        VG.resize(firstM);
        PsiR.resize(firstM, NRz);
        PsiS.resize(firstM, NSz);
        dPsidzR.resize(firstM, NRz);
        dPsidzS.resize(firstM, NSz);
    }

    void setZero()
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

// @brief 初始化结构体parameters
struct parameters
{
    // @brief 标题
    std::string Title;

    // @brief 任务数（声源*发射声线个数）
    int totalTasks;

    // 距离剖面个数
    int NProf;
    // @brief 距离剖面向量
    VectorXd RProf;
    int NMeshMax; // 最大网格数
    int NMediaMax;   // 最大媒质数

    // @brief 频率信息
    std::unique_ptr<FreqInfo> freqinfo;

    // // @brief 媒质数
    // int NMedia;

    // @brief 吸收单位
    Atten_Mode AttenUnit;

    // @brief 声源位置
    std::unique_ptr<Position> Pos;

    // @brief 声速剖面参数
    std::unique_ptr<SSPStructure []> SSP;

    // @brief 边界参数
    std::unique_ptr<BdryType> Bdry;

    // @brief 边界形状参数
    // Matrix<BdryPt, 1, Dynamic> Top, Bot;

    // @brief 顶部边界点数
    // int NatiPts;

    ReflectionCoefInfo ReflectionCoef; // @brief 反射系数参数

    // @brief 底部边界点数
    // int NbtyPts;

    // bool isTopSet, isBotSet; // @brief 是否设置了顶部和底部边界

    // @brief 实例化SrcBmPat为SBP
    std::unique_ptr<SrcBmPat> SBP;

    // @brief 半空间为'g'时吸声系数单位
    Atten_Mode BG_AttenUnit;

    int Number_to_Echo = 21;

    // @brief epsilon
    std::complex<double> epsilon;

    bool is_Velocity = false; // @brief 是否计算振速
    kkc_Log *log;             // 日志

    // 核心参数
    double cLow;            // 最小相速度
    double cHigh;           // 最大相速度
    double Rmax;            // 最大距离
    Source_Mode SourceType; // 声源模式

    // 网格参数
    MeshParams mesh; // 网格参数

    // 半空间参数
    std::unique_ptr<HSInfo []> HSTop; // 顶部半空间
    std::unique_ptr<HSInfo []> HSBot; // 底部半空间

    // 运行模式
    Run_Mode runMode; // Kraken运行模式

    // 相干和非相干
    CoherenceType coherenceType; // 相干和非相干

    // 模式类型
    ModeType modeType; // 模式类型

    // 输出控制
    bool outputModes; // 是否输出模式
    bool outputField; // 是否输出声场
};

// @brief 三对角矩阵结构体，kraken计算重要的中间变量
struct TridMtx
{
    // @brief 深度向量
    VectorXd z;
    // @brief 声速剖面插值点声速
    VectorXcd cp_int;
    VectorXcd cs_int;
    VectorXd rho_int;



    // 计算矩阵
    VectorXd B1;
    VectorXd B1C;
    VectorXd B2;
    VectorXd B3;
    VectorXd B4;
    VectorXd rho;

    // 本征函数
    VectorXd psi;     // 原始mesh的本征函数值
    VectorXd dpsidz;     // 原始mesh的本征函数值的z方向导数

    // @brief 网格参数
    VectorXi N;   // 各层的网格点数
    VectorXd h;   // 各层的网格步长
    VectorXd hV;  // 网格步长向量
    VectorXi Loc; // 各层的网格点索引

    // @brief 相速度范围
    double cLow;
    double cHigh;

    void resize(int Maxsize, int NMediaMax, int NSets)
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
struct kkc_output
{
    // 本征值参数
    EigenParams *eigen;                // 本征值参数
    std::complex<float> *u_AllSources; // 声压，一维化存储，内存优化较好
    std::complex<float> *v_AllSources; // 垂直振速，一维化存储，内存优化较好
    std::complex<float> *h_AllSources; // 水平振速，一维化存储，内存优化较好
};

#endif // PARAMS_H
