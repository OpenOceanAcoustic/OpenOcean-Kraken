#ifndef Kraken_JSON_IN_OUT_HPP
#define Kraken_JSON_IN_OUT_HPP

#include "OpenOceanKrakenParams.h"
#include "OpenOceanKrakenInterface.h"

namespace OpenOceanKraken
{   

    //各类enum
    void to_json(OpenOcean_json& out,const SSP_Mode& mode);
    void from_json(const OpenOcean_json& in,SSP_Mode& mode);

    void to_json(OpenOcean_json& out,const Media_Mode& mode);
    void from_json(const OpenOcean_json& in,Media_Mode& mode);

    void to_json(OpenOcean_json& out,const AttenuationUnit& mode);
    void from_json(const OpenOcean_json& in,AttenuationUnit& mode);

    void to_json(OpenOcean_json& out,const OceanAbsorptionModel& mode);
    void from_json(const OpenOcean_json& in,OceanAbsorptionModel& mode);

    void to_json(OpenOcean_json& out,const Atten_Mode& mode);
    void from_json(const OpenOcean_json& in,Atten_Mode& mode);

    void to_json(OpenOcean_json& out,const BC_Mode& mode);
    void from_json(const OpenOcean_json& in,BC_Mode& mode);

    void to_json(OpenOcean_json& out,const Source_Mode& mode);
    void from_json(const OpenOcean_json& in,Source_Mode& mode);

    void to_json(OpenOcean_json& out,const Grid_Mode& mode);
    void from_json(const OpenOcean_json& in,Grid_Mode& mode);

    void to_json(OpenOcean_json& out,const Run_Mode& mode);
    void from_json(const OpenOcean_json& in,Run_Mode& mode);

    void to_json(OpenOcean_json& out,const CoherenceType& mode);
    void from_json(const OpenOcean_json& in,CoherenceType& mode);

    void to_json(OpenOcean_json& out,const ModeType& mode);
    void from_json(const OpenOcean_json& in,ModeType& mode);

    //HSInfo
    void to_json(OpenOcean_json& out,const HSInfo& HS);
    void from_json(const OpenOcean_json& in,HSInfo& HS);

    //position
    void to_json(OpenOcean_json& out,const Position& pos);
    void from_json(const OpenOcean_json& in,Position& pos);

    //Reflection
    void to_json(OpenOcean_json& out,const ReflectionCoef& ref);
    void from_json(const OpenOcean_json& in,ReflectionCoef& ref);
    //ReflectionCoefInfo
    void to_json(OpenOcean_json& out,const ReflectionCoefInfo& ref);
    void from_json(const OpenOcean_json& in,ReflectionCoefInfo& ref);

    //SrcBmPat
    void to_json(OpenOcean_json& out,const SrcBmPat& pat);
    void from_json(const OpenOcean_json& in,SrcBmPat& pat);

    //freqInfo
    void to_json(OpenOcean_json& out,const FreqInfo& freq);
    void from_json(const OpenOcean_json& in,FreqInfo& freq);

    //params
    void to_json(OpenOcean_json& out,const OOK_parameters& params);
    void from_json(const OpenOcean_json& in,OOK_parameters& params);
    





    
    //ssp命名空间里面的转换
    namespace ssp
    {

        //SSPLayer
        void to_json(OpenOcean_json& out,const SSPLayer& Layer);
        void from_json(const OpenOcean_json& in,SSPLayer& Layer);

        //Range_Independent_Area
        void to_json(OpenOcean_json& out,const Range_Independent_Area& Area);
        void from_json(const OpenOcean_json& in,Range_Independent_Area& Area);
    }



}











#endif
