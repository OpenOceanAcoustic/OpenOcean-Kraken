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
    Media_Mode Material;
    double sigma;
    // @brief 声速剖面点数
    int NPts;
    // @brief 声速剖面细分点数
    int N;
    // 层厚度
    double depth;
    // 细分步长
    double h;
    // 介质层定义
    double beta;
    double ft;

    // @brief 深度向量
    VectorXd z;

    VectorXd alphaR; // 声速，纵波速度
    VectorXd alphaI; // 横波速度
    VectorXd betaR;  // 纵波衰减
    VectorXd betaI;  // 横波衰减
    // @brief 密度向量
    VectorXd rho;

    // @brief 声速向量
    VectorXcd cp;
    VectorXcd cs;
    VectorXcd cp_int;
    VectorXcd cs_int;
    VectorXd rho_int;

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

// kraken 计算矩阵
struct KrakenMatrix
{
    VectorXd B1;
    VectorXd B1C;
    VectorXd B2;
    VectorXd B3;
    VectorXd B4;
    VectorXd rho;
    int modeCount;
};


// 有限差分网格参数
struct MeshParams
{
    VectorXi Loc;
    int NSets = 5;   // 网格集数量
    int NV[5] = {1, 2, 4, 8, 16};   // Richardson外推系数数组
    VectorXi N;  // 各层的网格点数
    VectorXd h;  // 各层的网格步长
    VectorXd hV; // 网格步长向量
};

// 本征值相关参数
struct EigenParams
{
    int M;             // 模式数量
    VectorXd EVMat;    // 本征值矩阵(一维向量化)
    VectorXd Extrap;   // 外推矩阵
    VectorXcd k;        // 波数向量
    VectorXd VG;       // 群速度向量
    int LRecordLength; // 记录长度
    int IRecProfile;   // 记录指针
};

// 本征函数结构体
struct EigenFunction
{
    MatrixXcd phi;  // 本征函数值
    VectorXi modes; // 模式索引
    VectorXd depth; // 深度向量
};

// Kraken特有的运行模式
enum class Run_Mode
{
    MODE_M_Modes, // 只计算简正波模式
    MODE_F_Field, // 计算声场
    MODE_B_Both   // 同时计算模式和声场
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
    size_t totalTasks;

    // @brief 频率信息
    FreqInfo *freqinfo;

    // @brief 媒质数
    int NMedia;

    // @brief 声速剖面类型
    SSP_Mode SSPType;
    // @brief 吸收单位
    Atten_Mode AttenUnit;

    // @brief 声源位置
    Position *Pos;

    int FirstAcoustic;
    int LastAcoustic;

    // @brief 声速剖面参数
    SSPStructure *SSP;

    // @brief 边界参数
    BdryType *Bdry;

    // @brief 边界形状参数
    Matrix<BdryPt, 1, Dynamic> Top, Bot;

    // @brief 顶部边界点数
    int NatiPts;

    ReflectionCoefInfo ReflectionCoef; // @brief 反射系数参数

    // @brief 底部边界点数
    int NbtyPts;

    bool isTopSet, isBotSet; // @brief 是否设置了顶部和底部边界

    // @brief 实例化SrcBmPat为SBP
    SrcBmPat *SBP;

    // @brief 半空间为'g'时吸声系数单位
    Atten_Mode BG_AttenUnit;

    int Number_to_Echo = 21;

    // @brief epsilon
    std::complex<double> epsilon;

    // @brief iBeamWindow2
    double iBeamWindow2;

    // @brief RadMax
    double RadMax;

    // @brief ft
    double ft;

    bool is_Velocity = false; // @brief 是否计算振速
    kkc_Log *log;             // 日志

    // 核心参数
    double Clow;            // 最小相速度
    double Chigh;           // 最大相速度
    double Rmax;            // 最大距离
    Source_Mode SourceType; // 声源模式

    // 网格参数
    MeshParams mesh; // 网格参数

    // 半空间参数
    HSInfo HSTop; // 顶部半空间
    HSInfo HSBot; // 底部半空间

    // 运行模式
    Run_Mode runMode; // Kraken运行模式

    // 模式类型
    ModeType modeType; // 模式类型

    // 输出控制
    bool outputModes; // 是否输出模式
    bool outputField; // 是否输出声场
};

// 输出结构
struct kkc_output
{
    // 本征值参数
    EigenParams eigen;                 // 本征值参数
    std::complex<float> *u_AllSources; // 声压，一维化存储，内存优化较好
    std::complex<float> *v_AllSources; // 垂直振速，一维化存储，内存优化较好
    std::complex<float> *h_AllSources; // 水平振速，一维化存储，内存优化较好
};

#endif // PARAMS_H
