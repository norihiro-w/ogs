/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file PorousMedia.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include "BaseLib/CodingTools.h"


namespace NumLib
{
class ITXFunction;
}

namespace MaterialLib
{

struct PorousMedia
{
    NumLib::ITXFunction* hydraulic_conductivity;
    NumLib::ITXFunction* permeability;
    NumLib::ITXFunction* porosity;
    NumLib::ITXFunction* storage;
    NumLib::ITXFunction* geo_area;
    NumLib::ITXFunction* dispersivity_long;
    NumLib::ITXFunction* dispersivity_trans; 

    PorousMedia()
    {
        BaseLib::zeroObject(
                hydraulic_conductivity,
                permeability,
                porosity,
                storage,
                geo_area
                );
    }

    ~PorousMedia()
    {
        BaseLib::releaseObject(
                hydraulic_conductivity,
                permeability,
                porosity,
                storage,
                geo_area
                );
    }
};

} //end
 
