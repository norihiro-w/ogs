/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 */

#pragma once

#include <cassert>
#include <vector>

#include "MathLib/DataType.h"
#include "Models.h"

namespace MaterialLib
{

template <typename T>
class TemplateBasicModel
{
public:
    enum class Type
    {
        Constant,
        Linear
    };

    TemplateBasicModel(Type type, double parameter)
    :_model_type(type), _parameters(1, parameter)
    {
    }

    TemplateBasicModel(Type type, std::vector<double> &parameters)
    :_model_type(type), _parameters(parameters)
    {
    }

    T operator()(const StateVariables* variables)
    {
        switch(_model_type)
        {
        case Type::Constant:                   // rho = const
            {
                return _parameters[0];
            }
            break;
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
                return rho_0 *
                          (1.   + drho_dp * (variables->p - p_0)
                                + drho_dC * (variables->C - C_0)
                                + drho_dT * (variables->T - T_0));
            }
            break;
        }
        return 0;
    }

private:
    Type _model_type;
    std::vector<double> _parameters;
};

class BasicModelTensor
{
public:
    enum class Type
    {
        Constant
    };

    BasicModelTensor(Type type, unsigned dim, double parameter)
    :_model_type(type), _dim(dim), _parameters(1, parameter)
    {
    }

    BasicModelTensor(Type type, unsigned dim, std::vector<double> &parameters)
    :_model_type(type), _dim(dim), _parameters(parameters)
    {
    }

    void operator()(const StateVariables* /*variables*/, unsigned ele_dim, unsigned global_dim, const MathLib::RotationMatrix* matR, MathLib::LocalMatrix &ret)
    {
        _dim = ele_dim; //TODO
        switch(_model_type)
        {
        case Type::Constant:                   // rho = const
            {
                ret.setIdentity(_dim, _dim);
                ret *= _parameters[0];
            }
            break;
        default:
            break;
        }
        to_global(ret, global_dim, matR);
    }

private:
    void to_global(MathLib::LocalMatrix &local, unsigned global_dim, const MathLib::RotationMatrix* matR)
    {
        if (local.rows() < global_dim) {
            assert(matR!=nullptr);
            MathLib::LocalMatrix local2 = MathLib::LocalMatrix::Zero(global_dim, global_dim);
            local2.block(0, 0, local.rows(), local.cols()) = local.block(0, 0, local.rows(), local.cols());
            local.noalias() = (*matR) * local2 * matR->transpose();
        }
    }

    Type _model_type;
    unsigned _dim;
    std::vector<double> _parameters;
};

} //end
