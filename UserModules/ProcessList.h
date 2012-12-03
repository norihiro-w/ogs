/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file ProcessList.h
 *
 * Created on 2012-07-04 by Norihiro Watanabe
 */

#pragma once

//------------------------------------------------------------------------------
// List of active modules
// * This file should be called in body files (e.g. cpp)
//------------------------------------------------------------------------------
#include "FemGroundwaterFlow/Head.h"
#include "FemGroundwaterFlow/HeadToElementVelocity.h"
//#include "FemGroundwaterFlow/PressureToHead.h"
#include "FemGroundwaterFlow/PressureToElementVelocity.h"
//#include "FemGroundwaterFlow/LiquidPressure.h"
#include "FemMassTransport/Concentration.h"
//#include "FemHeatTransport/Temperature.h"
#include "FemKinReactGIA/Concentrations.h"
#include "FemDeformationTotalForm/Displacement.h"
//#include "FemDeformationTotalForm/ElementStressStrain.h"
//#include "FemDeformationTotalForm/NodalStressStrain.h"
//#include "FemPoroelastic/DisplacementPressure.h"
//#include "Xfem/XFEM_EXAMPLE_CRACK1.h"
//#include "THMmf/Pmf.h"
//#include "THMmf/Tmf.h"
//#include "THMmf/Umf.h"
//#include "ProcessBuilder.h"

//OGS_LINK_PROCESS(GROUNDWATER_FLOW, FunctionHead);
