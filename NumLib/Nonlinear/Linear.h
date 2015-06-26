/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file Linear.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include "INonlinearSolver.h"

namespace NumLib
{

/**
 * \brief Linear
 */
template <class F_LINEAR>
class Linear : public INonlinearSolver
{
    F_LINEAR* _linear_f;
public:
    explicit Linear(F_LINEAR* linear_f) : _linear_f(linear_f) {};
    virtual ~Linear() {};

    virtual void solve(MathLib::IVector &x)
    {
        (*_linear_f)(x);
    }
};

} //end

