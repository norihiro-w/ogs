/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file Picard.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include "MathLib/Nonlinear/Picard.h"
#include "INonlinearSolver.h"

namespace NumLib
{

/**
 * \brief Picard
 */
template <class F_LINEAR>
class Picard : public INonlinearSolver
{
public:
    Picard(F_LINEAR* linear_f) : _linear_f(linear_f)
    {
        _x_old = _dx = nullptr;
    };

    virtual ~Picard()
    {
        delete _x_old;
        delete _dx;
    };

    virtual void solve(MathLib::IVector &x)
    {
        if (_x_old == nullptr) {
            _x_old = x.duplicate();
            _dx = x.duplicate();
        }

        MathLib::Nonlinear::Picard picard;
        picard.setRelTolerance(getOption().error_tolerance);
        picard.setMaxIterations(getOption().max_iteration);
        picard.solve(*_linear_f, x);
    }

private:
    F_LINEAR* _linear_f;
    MathLib::IVector* _x_old;
    MathLib::IVector* _dx;
};

} //end

