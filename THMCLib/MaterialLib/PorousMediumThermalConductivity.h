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

class PorousMediumThermalConductivity
{
public:
    enum Type
    {
        Constant,
        ARITHMETIC
    };

    PorousMediumThermalConductivity(Type type, std::vector<double> &parameters)
    :_model_type(type), _parameters(parameters)
    {
    }

    void operator()(const StateVariables* variables, const PorousMediumProperty &pm, const SolidProperty &s, const std::vector<FluidProperty> &f, unsigned global_dim, MathLib::LocalMatrix &val)
    {
        switch(_model_type)
        {
        case Type::Constant:                   // rho = const
        {
            double rho_0 = _parameters[0];
            val.setIdentity(global_dim, global_dim);
            val *= rho_0;
        }
            break;

        case Type::ARITHMETIC:
        {
            val.noalias() = (1-pm.n)*s.lambda + pm.n*f[0].lambda*MathLib::LocalMatrix::Identity(global_dim, global_dim);
            break;
        }
        default:
            std::cout << "Error in PorousMediumThermalConductivity() no valid model" <<
            std::endl;
            break;
        }
    }

private:
    Type _model_type;
    std::vector<double> _parameters;
};

} //end
