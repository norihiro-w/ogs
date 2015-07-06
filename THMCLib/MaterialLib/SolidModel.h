/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 */

#pragma once


#include "BasicModel.h"
#include "MaterialProperties.h"

namespace MaterialLib
{

typedef TemplateBasicModel<double> SolidDensity;
typedef TemplateBasicModel<double> SolidThermalExpansion;
typedef TemplateBasicModel<double> PoissonRatio;
typedef TemplateBasicModel<double> YoungsModulus;
typedef TemplateBasicModel<double> SolidSpecificHeat;
typedef BasicModelTensor SolidThermalConductivity;

struct SolidModel
{
    SolidDensity* density;
    SolidThermalExpansion* thermal_expansion;
    PoissonRatio* poisson_ratio;
    YoungsModulus* Youngs_modulus;
    SolidSpecificHeat* specific_heat;
    SolidThermalConductivity* thermal_conductivity;


    SolidModel()
    : density(nullptr), thermal_expansion(nullptr), poisson_ratio(nullptr),
      Youngs_modulus(nullptr), specific_heat(nullptr), thermal_conductivity(nullptr)
    {
    }

    ~SolidModel()
    {
        delete density;
        delete thermal_expansion;
        delete poisson_ratio;
        delete Youngs_modulus;
        delete specific_heat;
        delete thermal_conductivity;
    }

    SolidProperty operator()(const StateVariables &var, const unsigned dim) const
    {
        SolidProperty v;
        if (density) v.rho = (*density)(&var);
        if (specific_heat) v.cp = (*specific_heat)(&var);
        if (thermal_conductivity) v.lambda = (*thermal_conductivity)(&var, dim);
        return v;
    }
};

}
