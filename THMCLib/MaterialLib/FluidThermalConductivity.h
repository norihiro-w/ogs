/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 */

#pragma once

#include "MathLib/DataType.h"
#include "Models.h"

namespace MaterialLib
{

class FluidThermalConductivity
{
public:
    enum Type
    {
        Constant,
        Linear,
        Curve,
        IAPWS_IF97
    };

    FluidThermalConductivity(Type type, std::vector<double> &parameters)
    :_model_type(type), _parameters(parameters)
    {
    }

    double operator()(const StateVariables* variables)
    {
        double density = 0;

        //----------------------------------------------------------------------
        switch(_model_type)
        {
        case Type::Curve:                   // rho = f(x)
//            int fct_number = 0;
//            int gueltig;
            //density = GetCurveValue(fct_number,0,variables->p,&gueltig);
            break;
        case Type::Constant:                   // rho = const
        {
            double rho_0 = _parameters[0];
            density = rho_0;
        }
            break;

        case Type::IAPWS_IF97:
        {
            density = SalineWaterThermalConductivity_IAPWS_IF97(variables->p,variables->T,variables->C);
            break;
        }
        case Type::Linear:
            {
                // rho(p,C,T) = rho_0*(1+beta_p*(p-p_0)+beta_C*(C-C_0)+beta_T*(T-T_0))
                unsigned i = 0;
                double rho_0 = _parameters[i++];
                double p_0 = _parameters[i++];
                double drho_dp = _parameters[i++];
                double T_0 = _parameters[i++];
                double drho_dT = _parameters[i++];
                double C_0 = _parameters[i++];
                double drho_dC = _parameters[i++];
                density = rho_0 *
                          (1.   + drho_dp * (variables->p - p_0)
                                + drho_dC * (variables->C - C_0)
                                + drho_dT * (variables->T - T_0));
            }
            break;
        default:
            std::cout << "Error in CFluidProperties::Density: no valid model" <<
            std::endl;
            break;
        }

        return density;
    }

private:
    Type _model_type;
    std::vector<double> _parameters;
};

} //end
