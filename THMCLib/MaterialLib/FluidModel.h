/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 */

#pragma once

#include "MathLib/DataType.h"
#include "Models.h"
#include "MaterialProperties.h"
#include "FluidDensity.h"
#include "FluidViscosity.h"
#include "FluidSpecificHeat.h"
#include "FluidThermalConductivity.h"

namespace MaterialLib
{

struct FluidModel
{
    FluidDensity* density;
    FluidViscosity* viscosity;
    FluidSpecificHeat* specific_heat;
    FluidThermalConductivity* thermal_conductivity;

    FluidModel() : density(nullptr), viscosity(nullptr), specific_heat(nullptr), thermal_conductivity(nullptr)
    {}

    ~FluidModel()
    {
        delete density;
        delete viscosity;
        delete specific_heat;
        delete thermal_conductivity;
    }

    FluidProperty operator()(const StateVariables &var) const
    {
        FluidProperty v;
        if (density) v.rho = (*density)(&var);
        if (viscosity) v.mu = (*viscosity)(&var);
        if (specific_heat) v.cp = (*specific_heat)(&var);
        if (thermal_conductivity) v.lambda = (*thermal_conductivity)(&var);
        return v;
    }
};

} //end
