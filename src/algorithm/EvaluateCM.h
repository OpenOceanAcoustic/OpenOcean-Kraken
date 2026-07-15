#ifndef EVALUATE_CM_H
#define EVALUATE_CM_H

#include "OpenOceanKrakenParams.h"

namespace OpenOceanKraken
{
    void EvaluateCM(const EigenParams *profiles,
                    int profileCount,
                    const OOK_parameters &params,
                    int isz,
                    std::complex<float> *uAllSources);
}

#endif
