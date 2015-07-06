/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 */

#pragma once

#include "MathLib/DataType.h"

namespace MaterialLib
{

struct StateVariables
{
    double p;
    double T;
    double C;
};

double SalineWaterDensityIAPWS_IF97(double Press, double TempK, double Conc);
double SalineWaterSpecificHeat_IAPWS_IF97(double Press, double TempK, double Conc);
double SalineWaterThermalConductivity_IAPWS_IF97(double Press, double TempK, double Conc);

double LiquidViscosity_HP(double Press,double TempK,double C_0);

/**
 * Dynamic viscosity of high-concentration salt water based on Lever&Jackson(1985),
 * Hassanizadeh(1988), and Mercer&Pinder(1974)
 *
 * Reference: FEFLOW Reference manual pg 31, (1-100)
 */
double LiquidViscosity_LJH_MP1(double c,double T, double rho, double my_0);

/**
 * Dynamic viscosity of high-concentration salt water based on Lever&Jackson(1985),
 * Hassanizadeh(1988), and Mercer&Pinder(1974)
 *
 * Reference: FEFLOW Reference manual pg 31, (1-100)
 */
double LiquidViscosity_LJH_MP2(double c,double T, double rho, double rho0, double my_C0, double my_T0, double my_0);


} //end
