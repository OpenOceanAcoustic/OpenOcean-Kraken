#ifndef EVALUATE_AD_H
#define EVALUATE_AD_H

#include "OpenOceanKrakenParams.h"

namespace OpenOceanKraken
{
    void EvaluateAD(const EigenParams *profiles,
                    int profileCount,
                    const OOK_parameters &params,
                    int isz,
                    std::complex<float> *uAllSources);
}

#endif
