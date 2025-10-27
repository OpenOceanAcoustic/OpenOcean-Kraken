#ifndef PCHIPMOD_H
#define PCHIPMOD_H


#include <cmath>

#include "splinec.h"


// 计算PCHIP的成员函数
void PCHIP(VectorXd& x, VectorXcd& y, int N, MatrixXcd &PolyCoef,MatrixXcd &csWork);

void h_del(VectorXd& x, VectorXcd& y, int ix, double &h1, double &h2, std::complex<double> &del1, std::complex<double> &del2);

std::complex<double> fprime_interior_Cmplx(const std::complex<double> &del1, const std::complex<double> &del2, const std::complex<double> &fprime);

std::complex<double> fprime_left_end_Cmplx(const std::complex<double> &del1, const std::complex<double> &del2, const std::complex<double> &fprime);


std::complex<double> fprime_right_end_Cmplx(const std::complex<double> &del1, const std::complex<double> &del2, const std::complex<double> &fprime);


double fprime_interior(double del1, double del2, double fprime);


double fprime_left_end(double del1, double del2, double fprime);


double fprime_right_end(double del1, double del2, double fprime);



#endif // PCHIPMOD_H