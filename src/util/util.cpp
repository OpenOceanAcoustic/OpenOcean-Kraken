#include "util.h"

// 线性插值
void linspace(double a, double b, int n, Eigen::VectorXd &x)
{
    // if (n < 2) {
    //     std::cout<<"linspace: Number of elements must be at least 2"<<std::endl;
    //     return;
    // }
    if (n < 1)
    {
        std::cout << "linspace: Number of elements must be at least 1" << std::endl;
        return;
    }

    x.resize(n);
    if (n == 1)
    {
        x[0] = a;
        return;
    }

    double h = (b - a) / (n - 1);
    for (int i = 0; i < n - 1; ++i)
    {
        x[i] = a + i * h;
    }
    x[n - 1] = b; // 手动设置最后一个元素为 b
}

// 插值(废弃)
void SubTab(VectorXd &x, int Nx)
{
    double deltax;

    if (Nx >= 3)
    {
        if (x[2] == -999.9)
        {
            if (x[1] == -999.9)
            {
                x[1] = x[0];
            }
            deltax = (x[1] - x[0]) / (Nx - 1);
            for (int i = 0; i < Nx; ++i)
            {
                x[i] = x[0] + i * deltax;
            }
        }
    }
}

// 某种排序
void Sort(VectorXd &x, int N)
{
    if (N == 1)
        return;

    double xTemp;
    int IRight, ILeft, IMiddle;

    for (int I = 1; I < N; ++I)
    {
        xTemp = x[I];

        if (xTemp < x[0])
        {
            for (int j = I; j > 0; --j)
            {
                x[j] = x[j - 1];
            }
            x[0] = xTemp;
        }
        else if (xTemp < x[I - 1])
        {
            IRight = I - 1;
            ILeft = 0;

            while (IRight > ILeft + 1)
            {
                IMiddle = (ILeft + IRight) / 2;
                if (xTemp < x[IMiddle])
                {
                    IRight = IMiddle;
                }
                else
                {
                    ILeft = IMiddle;
                }
            }

            for (int j = I; j > IRight + 1; --j)
            {
                x[j] = x[j - 1];
            }
            x[IRight] = xTemp;
        }
    }
}

double spacing(double const &x)
{
    return std::abs(std::nextafter(x, 0.0) - x);
}

bool isSmallValue(double delta, double ref, bool const_thresh)
{
    const double threadshold = const_thresh ? float(1.0e-4) : double(1.0e3) * spacing(ref);
    return delta < threadshold;
}

std::complex<float> cpxd2cpxf(const std::complex<double> &a)
{
    return std::complex<float>(float(a.real()), float(a.imag()));
}

std::complex<double> cpxf2cpxd(const std::complex<float> &a)
{
    return std::complex<double>(double(a.real()), double(a.imag()));
}
