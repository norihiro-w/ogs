/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 */

#pragma once

#include <vector>

#include "MathLib/DataType.h"
#include "Models.h"
#include "FluidDensity.h"

namespace MaterialLib
{

class FluidViscosity
{
public:
    enum Type
    {
        Constant,
        Linear,
        Curve,
        HP,
        LJH_MP1,
        LJH_MP2
    };

    FluidViscosity(Type type, std::vector<double> &parameters)
    :_viscosity_model(type), _parameters(parameters)
    {
    }

    double operator()(const StateVariables* variables)
    {
        double viscosity = 0.0;

        switch(_viscosity_model)
        {
        case Type::Curve:
//            int fct_number = 0;
//            int gueltig;
//            viscosity = GetCurveValue(fct_number,0,variables->p,&gueltig);
            break;
        case Type::Constant:                               // my = const
            viscosity = _parameters[0];
            break;
        case Type::Linear:                               // my(p) = my_0*(1+gamma_p*(p-p_0))
        {
            double my_0 = _parameters[0];
            double p_0 = _parameters[0];
            double dmy_dp = _parameters[0];
            viscosity = my_0 * (1. + dmy_dp * (variables->p - p_0));
        }
            break;
        case Type::HP:                               // my(p,C,T),
            viscosity = LiquidViscosity_HP(variables->p, variables->T, variables->C);
            break;
        case Type::LJH_MP1:
            {
                viscosity = LiquidViscosity_LJH_MP1(variables->C, variables->T, 0.0, 0.0); //c[g/L], T[C]
            }
            break;
        case Type::LJH_MP2:
        {
            FluidDensity* density = nullptr; //TODO
            double rho = (*density)(variables);
            double rho0=0, my_C0=0, my_T0=0, my_0=0;
            viscosity = LiquidViscosity_LJH_MP2(variables->C, variables->T, rho, rho0, my_C0, my_T0, my_0); //c[g/L], T[C]
        }
            break;
        }

        return viscosity;
    }

private:
    Type _viscosity_model;
    std::vector<double> _parameters;
};

} //end
