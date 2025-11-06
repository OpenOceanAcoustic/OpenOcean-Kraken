//沉积层参数计算模块
#ifndef ATTENMOD_H
#define ATTENMOD_H
#include "kkc_params.h"

#include <cmath>
std::complex<double> CRCI(double& z, double& c, double& alpha, double& freq, double& freq0, 
    const Atten_Mode& AttenUnit, double& beta, double& ft);
double Franc_Garr(double f);

bool isThorpe(const Atten_Mode& AttenUnit); //判断是否为Thorpe沉积层模型
bool isFranc_Garr(const Atten_Mode& AttenUnit);//判断是否为Franc-Garr沉积层模型

bool is_F(const Atten_Mode& AttenUnit);
bool is_L(const Atten_Mode& AttenUnit);
bool is_M(const Atten_Mode& AttenUnit);
bool is_m(const Atten_Mode& AttenUnit);
bool is_N(const Atten_Mode& AttenUnit);
bool is_Q(const Atten_Mode& AttenUnit);
bool is_W(const Atten_Mode& AttenUnit);


#endif