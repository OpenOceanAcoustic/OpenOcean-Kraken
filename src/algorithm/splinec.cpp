#include "splinec.h"

void CSpline(VectorXd &TAU, MatrixXcd &C, int N, int IBCBEG, int IBCEND, int NDIM)
{
    int L = N - 1;
    std::complex<double> G, DTAU, DIVDF1, DIVDF3;

    for (int M = 1; M < N; ++M)
    {
        C(2, M) = TAU(M) - TAU(M-1);
        C(3, M) = (C(0, M) - C(0, M-1)) / C(2, M);
    }

    if (IBCBEG == 0)
    {
        if (N > 2)
        {
            C(3,0) = C(2,2);
            C(2,0) = C(2,1) + C(2,2);
            C(1,0) = ((C(2,1) + 2.0 * C(2,0)) * C(3,1) * C(2,2) + SQ(C(2,1)) * C(3,2)) / C(2,0);
        }
        else
        {
            C(3,0) = std::complex<double>(1.0, 0.0);
            C(2,0) = std::complex<double>(1.0, 0.0);
            C(1,0) = 2.0 * C(3,1);
        }
    }
    else if (IBCBEG == 1)
    {
        C(3,0) = std::complex<double>(1.0, 0.0);
        C(2,0) = std::complex<double>(0.0, 0.0);
    }
    else if (IBCBEG == 2)
    {
        C(3,0) = std::complex<double>(2.0, 0.0);
        C(2,0) = std::complex<double>(1.0, 0.0);
        C(1,0) = 3.0 * C(3,1) - C(1,0) * C(2,1) / 2.0;
    }

    for (int M = 1; M < L; ++M)
    {
        G = -C(2,M+1) / C(3,M-1);
        C(1,M) = G * C(1,M-1) + 3.0 * (C(2,M) * C(3,M+1) + C(2,M+1) * C(3,M));
        C(3,M) = G * C(2,M-1) + 2.0 * (C(2,M) + C(2,M+1));
    }

    if (IBCEND != 1)
    {
        if (IBCEND == 0)
        {
            if ((N == 2 && IBCBEG == 0))
            {
                C(1,N-1) = C(3,N-1);
            }
            else if ((N == 3 && IBCBEG == 0) || N == 2)
            {
                C(1,N-1) = 2.0 * C(3,N-1);
                C(3,N-1) = complex<double>(1.0, 0.0);
                G = -1.0 / C(3,N-2);
            }
            else
            {
                G = C(2,N-2) + C(2,N-1);
                C(1,N-1) = ((C(2,N-1) + 2.0 * G) * C(3,N-1) * C(2,N-2) +
                               SQ(C(2,N-1)) * (C(0,N-2) - C(0,N-3)) /
                                   C(2,N-2)) /
                              G;
                G = -G / C(3,N-2);
                C(3,N-1) = C(2,N-2);
            }
        }
        else if (IBCEND == 2)
        {
            C(1,N-1) = 3.0 * C(3,N-1) + C(1,N-1) * C(2,N-1) / 2.0;
            C(3,N-1) = 2.0;
            G = -1.0 / C(3,N-2);
        }

        if (IBCBEG > 0 || N > 2)
        {
            C(3,N-1) = G * C(2,N-2) + C(3,N-1);
            C(1,N-1) = (G * C(1,N-2) + C(1,N-1)) / C(3,N-1);
        }
    }

    for (int J = L - 1; J >= 0; --J)
    {
        C(1,J) = (C(1,J) - C(2,J) * C(1,J+1)) / C(3,J);
    }

    for (int I = 1; I < N; ++I)
    {
        DTAU = C(2,I);
        DIVDF1 = (C(0,I) - C(0,I-1)) / DTAU;
        DIVDF3 = C(1,I-1) + C(1,I) - 2.0 * DIVDF1;
        C(2,I-1) = 2.0 * (DIVDF1 - C(1,I-1) - DIVDF3) / DTAU;
        C(3,I-1) = DIVDF3 / DTAU * (6.0 / DTAU);
    }

    C(2,N-1) = C(2,L-1) + (TAU(N-1) - TAU(L-1)) * C(3,L-1);

    C(3,N-1) = std::complex<double>(0.0, 0.0);
    for (int I = 0; I < L; ++I)
    {
        DTAU = TAU(I+1) - TAU(I);
        C(3,N-1) += DTAU * (C(0,I) + DTAU * (C(1,I) / 2.0 +
                                                 DTAU * (C(2,I) / 6.0 + DTAU * C(3,I) / 24.0)));
    }
    C(3,N-1) /= (TAU(N-1) - TAU(0));
}

void VSpline(VectorXd &TAU, VectorXcd &C, int M, int MDIM, VectorXcd &F, int N)
{
    int J = 0; // Start J at 0 instead of 1

    for (int I = 0; I < N; ++I)
    {
        int J1 = J + 1;
        // Check if the current point is within the next interval
        while (TAU[J1] < real(F[I]) && J1 < M)
        {
            ++J;
            J1 = J + 1;
        }
        // Calculate the distance from the start of the interval
        double H = real(F[I]) - TAU[J];
        // Use SPLINE function to interpolate
        F[I] = spline(&C[J * 4], H);
    }
}

std::complex<double> spline(const std::complex<double> *C, double H)
{
    std::complex<double> SPLINE;
    SPLINE = C[0] + H * (C[1] + H * (C[2] / 2.0 + H * C[3] / 6.0));
    return SPLINE;
}

// std::complex<double> splinex(const std::complex<double> *C, double H)
// {
//     std::complex<double> SPLINEX;
//     SPLINEX = C[1] + H * (C[2] + H * C[3] / 2.0);
//     return SPLINEX;
// }

// std::complex<double> splinexx(const std::complex<double> *C, double H)
// {
//     std::complex<double> SPLINEXX;
//     SPLINEXX = C[2] + H * C[3];
//     return SPLINEXX;
// }

void SplineALL(MatrixXcd &C, int &iSegz, double &H, std::complex<double> &F, std::complex<double> &FX, std::complex<double> &FXX){
    const double HALF = 0.5;
    const double SIXTH = 1.0 / 6.0;

    F=C(0,iSegz) +H*(C(1,iSegz) +H*(HALF*C(2,iSegz) +SIXTH*H*C(3,iSegz)));
    FX=C(1,iSegz) +H*(C(2,iSegz) +H*HALF*C(3,iSegz));
    FXX=C(2,iSegz) +H*C(3,iSegz);

}