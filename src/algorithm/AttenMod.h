// 沉积层参数计算模块
#ifndef ATTENMOD_H
#define ATTENMOD_H
#include "OpenOceanKrakenParams.h"

#include <cmath>

namespace OpenOceanKraken
{   

    // @brief 在基础衰减上叠加海洋吸收模型
    double addOceanAbsorption(double alphaT, double freq, const Atten_Mode& AttenUnit);
    double parseAttenuation(double freq, double freq0,double alpha, double ft, double beta, double c,const Atten_Mode& AttenUnit);

    std::complex<double> CRCI(double &z, double &c, double &alpha, double &freq, double &freq0,
                              const Atten_Mode &AttenUnit, double &beta, double &ft);

    double Franc_Garr(double f);

    
    bool isThorpe(const Atten_Mode &AttenUnit);     // 判断是否为Thorpe沉积层模型
    bool isFranc_Garr(const Atten_Mode &AttenUnit); // 判断是否为Franc-Garr沉积层模型
    bool is_F(const Atten_Mode &AttenUnit);
    bool is_L(const Atten_Mode &AttenUnit);
    bool is_M(const Atten_Mode &AttenUnit);
    bool is_m(const Atten_Mode &AttenUnit);
    bool is_N(const Atten_Mode &AttenUnit);
    bool is_Q(const Atten_Mode &AttenUnit);
    bool is_W(const Atten_Mode &AttenUnit);

}

#endif