#include "RefCoef.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace OpenOceanKraken
{
    namespace
    {
        std::complex<double> PolyC(const std::complex<double> &x0,
                                   const std::vector<std::complex<double>> &x,
                                   const std::vector<std::complex<double>> &f)
        {
            const size_t n = f.size();
            std::vector<std::complex<double>> ft = f;
            std::vector<std::complex<double>> h(n);
            for (size_t i = 0; i < n; ++i)
            {
                h[i] = x[i] - x0;
            }

            for (size_t i = 1; i < n; ++i)
            {
                for (size_t j = 0; j + i < n; ++j)
                {
                    ft[j] = ft[j] + h[j] * (ft[j] - ft[j + 1]) / (h[j + i] - h[j]);
                }
            }
            return ft[0];
        }
    }

    void InterpolateReflectionCoefficient(ReflectionCoef &RInt, const Eigen::Matrix<ReflectionCoef, 1, Eigen::Dynamic> &R)
    {

        int iLeft = 1;
        int iRight = R.size();
        double thetaIntr = RInt.theta;

        if (thetaIntr < R(iLeft - 1).theta)
        {
            RInt.R = 0.0;
            RInt.phi = 0.0;
        }
        else if (thetaIntr > R(iRight - 1).theta)
        {
            RInt.R = 0.0;
            RInt.phi = 0.0;
        }
        else
        {
            while (iLeft != iRight - 1)
            {
                int iMid = (iLeft + iRight) / 2;
                if (R(iMid - 1).theta > thetaIntr)
                {
                    iRight = iMid;
                }
                else
                {
                    iLeft = iMid;
                }
            }

            double alpha = (RInt.theta - R(iLeft - 1).theta) / (R(iRight - 1).theta - R(iLeft - 1).theta);
            RInt.R = (1 - alpha) * R(iLeft - 1).R + alpha * R(iRight - 1).R;
            RInt.phi = (1 - alpha) * R(iLeft - 1).phi + alpha * R(iRight - 1).phi;
        }
    }

    void InterpolateIRC(const std::complex<double> &x, std::complex<double> &f, std::complex<double> &g,
                        int &iPower, const InternalReflectionCoefInfo &irc)
    {
        const int nk = static_cast<int>(irc.xTab.size());
        if (!irc.isSet || nk == 0 || irc.fTab.size() != nk || irc.gTab.size() != nk || irc.iTab.size() != nk)
        {
            throw std::runtime_error("Internal reflection coefficient table is not available.");
        }

        const double xReal = std::real(x);
        if (nk == 1 || xReal < irc.xTab(0))
        {
            f = irc.fTab(0);
            g = irc.gTab(0);
            iPower = irc.iTab(0);
            return;
        }
        if (xReal > irc.xTab(nk - 1))
        {
            f = irc.fTab(nk - 1);
            g = irc.gTab(nk - 1);
            iPower = irc.iTab(nk - 1);
            return;
        }

        int iLeft = 0;
        int iRight = nk - 1;
        while (iLeft != iRight - 1)
        {
            const int iMid = (iLeft + iRight) / 2;
            if (irc.xTab(iMid) > xReal)
            {
                iRight = iMid;
            }
            else
            {
                iLeft = iMid;
            }
        }

        iLeft = std::max(iLeft - (3 - 2) / 2, 0);
        iRight = std::min(iRight + (3 - 1) / 2, nk - 1);

        const int nAct = iRight - iLeft + 1;
        std::vector<std::complex<double>> xT(nAct);
        std::vector<std::complex<double>> fT(nAct);
        std::vector<std::complex<double>> gT(nAct);
        for (int i = 0; i < nAct; ++i)
        {
            const int j = i + iLeft;
            const int iDel = irc.iTab(j) - irc.iTab(iLeft);
            const double scale = std::pow(10.0, iDel);
            xT[i] = std::complex<double>(irc.xTab(j), 0.0);
            fT[i] = irc.fTab(j) * scale;
            gT[i] = irc.gTab(j) * scale;
        }

        f = PolyC(x, xT, fT);
        g = PolyC(x, xT, gT);
        iPower = irc.iTab(iLeft);
    }

}
