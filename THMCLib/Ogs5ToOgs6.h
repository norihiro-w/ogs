/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file Ogs5ToOgs6.h
 *
 * Created on 2012-07-13 by Norihiro Watanabe
 */

#pragma once

#include <boost/property_tree/ptree.hpp>
#include <MaterialLib/SolidModel.h>

#include "FileIO/Legacy/Ogs5FemIO.h"
#include "MaterialLib/PorousMediumModel.h"
#include "GeoLib/GEOObjects.h"
#include "Ogs6FemData.h"

namespace ogs6
{

namespace Ogs5ToOgs6
{

void convertSolidProperty(const ogs5::CSolidProperties &msp, MaterialLib::SolidModel &solid);

void convertPorousMediumProperty(const ogs5::CMediumProperties &mmp, MaterialLib::PorousMediumModel &pm);

bool convert(const ogs5::Ogs5FemData &ogs5fem, THMCLib::Ogs6FemData &ogs6fem, boost::property_tree::ptree &option);
};
} //end
