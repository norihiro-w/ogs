/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 */

#pragma once


#include "BasicModel.h"

namespace MaterialLib
{

struct FluidProperty
{
    double rho;
    double mu;
    double cp;
    double lambda;
};

struct SolidProperty
{
    double rho;
    double mu;
    double cp;
    MathLib::LocalMatrix lambda;
};


struct PorousMediumProperty
{
    MathLib::LocalMatrix K; /// hydraulic conductivity
    MathLib::LocalMatrix k; /// permeability
    double n; /// porosity
    double Ss; /// specific storage
    double geo_area = 1.0; /// volume multiplier
    double C_alpha_L; /// solute dispersion long
    double C_alpha_T; /// solute dispersion trans
    double Cp; /// heat capacity
    MathLib::LocalMatrix lamba; /// heat conductivity
};
}
