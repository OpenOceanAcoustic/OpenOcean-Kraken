#include "pchipMod.h"

namespace OpenOceanKraken
{
    void PCHIP(Eigen::VectorXd &x, Eigen::VectorXcd &y, int N, Eigen::MatrixXcd &PolyCoef, Eigen::MatrixXcd &csWork)
    {
        int ix, iBCBeg, iBCEnd;
        double h, h1, h2;
        std::complex<double> del1, del2, f1, f2, f1prime, f2prime, fprimeT;

        if (N == 2)
        {
            PolyCoef(0, 0) = y(0);
            PolyCoef(1, 0) = (y(1) - y(0)) / (x(1) - x(0));
            PolyCoef(2, 0) = std::complex<double>(0.0, 0.0);
            PolyCoef(3, 0) = std::complex<double>(0.0, 0.0);
        }
        else
        {
            for (ix = 0; ix < N; ++ix)
            {
                PolyCoef(0, ix) = y(ix);
            }
            h_del(x, y, 1, h1, h2, del1, del2);
            fprimeT = ((2.0 * h1 + h2) * del1 - h1 * del2) / (h1 + h2);
            PolyCoef(1, 0) = fprime_left_end_Cmplx(del1, del2, fprimeT);
            h_del(x, y, N - 2, h1, h2, del1, del2);
            fprimeT = (-h2 * del1 + (h1 + 2.0 * h2) * del2) / (h1 + h2);
            PolyCoef(1, N - 1) = fprime_right_end_Cmplx(del1, del2, fprimeT);
            iBCBeg = 1;
            iBCEnd = 1;

            csWork.row(0) = PolyCoef.row(0);
            csWork(1, 0) = PolyCoef(1, 0);
            csWork(1, N - 1) = PolyCoef(1, N - 1);
            CSpline(x, csWork, N, iBCBeg, iBCEnd, N);

            for (ix = 1; ix < N - 1; ++ix)
            {
                h_del(x, y, ix, h1, h2, del1, del2);
                PolyCoef(1, ix) = fprime_interior_Cmplx(del1, del2, csWork(1, ix));
            }
            for (ix = 0; ix < N - 1; ++ix)
            {
                h = x(ix + 1) - x(ix);
                f1 = PolyCoef(0, ix);
                f2 = PolyCoef(0, ix + 1);
                f1prime = PolyCoef(1, ix);
                f2prime = PolyCoef(1, ix + 1);
                PolyCoef(2, ix) = (3.0 * (f2 - f1) - h * (2.0 * f1prime + f2prime)) / (h * h);
                PolyCoef(3, ix) = (h * (f1prime + f2prime) - 2.0 * (f2 - f1)) / (h * h * h);
            }
        }
    }

    void h_del(Eigen::VectorXd &x, Eigen::VectorXcd &y, int ix, double &h1, double &h2, std::complex<double> &del1, std::complex<double> &del2)
    {
        h1 = x(ix) - x(ix - 1);
        h2 = x(ix + 1) - x(ix);
        del1 = (y(ix) - y(ix - 1)) / h1;
        del2 = (y(ix + 1) - y(ix)) / h2;
    }

    std::complex<double> fprime_interior_Cmplx(const std::complex<double> &del1, const std::complex<double> &del2, const std::complex<double> &fprime)
    {
        double fprime_r = fprime_interior(std::real(del1), std::real(del2), std::real(fprime));
        double fprime_i = fprime_interior(std::imag(del1), std::imag(del2), std::imag(fprime));

        return std::complex<double>(fprime_r, fprime_i);
    }

    std::complex<double> fprime_left_end_Cmplx(const std::complex<double> &del1, const std::complex<double> &del2, const std::complex<double> &fprime)
    {
        double fprime_r = fprime_left_end(std::real(del1), std::real(del2), std::real(fprime));
        double fprime_i = fprime_left_end(std::imag(del1), std::imag(del2), std::imag(fprime));

        return std::complex<double>(fprime_r, fprime_i);
    }

    std::complex<double> fprime_right_end_Cmplx(const std::complex<double> &del1, const std::complex<double> &del2, const std::complex<double> &fprime)
    {
        double fprime_r = fprime_right_end(std::real(del1), std::real(del2), std::real(fprime));
        double fprime_i = fprime_right_end(std::imag(del1), std::imag(del2), std::imag(fprime));

        return std::complex<double>(fprime_r, fprime_i);
    }

    double fprime_interior(double del1, double del2, double fprime)
    {
        double fprime_interior;

        // Check if derivative is within the trust region, project into it if not
        if (del1 * del2 > 0.0)
        {
            // Adjacent secant slopes have the same sign, enforce monotonicity
            if (del1 > 0.0)
            {
                fprime_interior = std::min(std::max(fprime, 0.0), 3.0 * std::min(del1, del2));
            }
            else
            {
                fprime_interior = std::max(std::min(fprime, 0.0), 3.0 * std::max(del1, del2));
            }
        }
        else
        {
            // Force the interpolant to have an extremum here
            fprime_interior = 0.0;
        }

        return fprime_interior;
    }

    double fprime_left_end(double del1, double del2, double fprime)
    {
        double fprime_left_end = fprime;

        if (del1 * fprime <= 0.0)
        {
            // Set derivative to zero if the sign differs from the sign of secant slope
            fprime_left_end = 0.0;
        }
        else if ((del1 * del2 <= 0.0) && (std::abs(fprime) > std::abs(3.0 * del1)))
        {
            // Adjust derivative value to enforce monotonicity
            fprime_left_end = 3.0 * del1;
        }

        return fprime_left_end;
    }

    double fprime_right_end(double del1, double del2, double fprime)
    {
        double fprime_right_end = fprime;

        if (del2 * fprime <= 0.0)
        {
            // Set derivative to zero if the sign differs from the sign of secant slope
            fprime_right_end = 0.0;
        }
        else if ((del1 * del2 <= 0.0) && (std::abs(fprime) > std::abs(3.0 * del2)))
        {
            // Adjust derivative value to enforce monotonicity
            fprime_right_end = 3.0 * del2;
        }

        return fprime_right_end;
    }

}
