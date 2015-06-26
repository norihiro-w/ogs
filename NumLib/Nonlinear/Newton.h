/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file Newton.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include "BaseLib/CodingTools.h"
#include "MathLib/Nonlinear/NewtonRaphson.h"
//#include "MathLib/Nonlinear/NRIterationStepInitializerDummy.h"

#include "INonlinearSolver.h"

namespace NumLib
{

/**
 * \brief Newton-Raphson solver for discrete systems
 */
template <
    class F_R, 
    class F_J
//    class T_STEP_INIT=MathLib::NRIterationStepInitializerDummy
>
class NewtonRaphson : public INonlinearSolver
{
public:

    NewtonRaphson(F_R* f_r, F_J* f_J, MathLib::IMatrix* J, MathLib::IVector* r)
    : _f_r(f_r), _f_J(f_J), _J(J), _r(r) //, _nl_step_init(NULL)
    {
    }

    virtual ~NewtonRaphson()
    {}

    virtual void solve(MathLib::IVector &x)
    {
        MathLib::Nonlinear::NewtonRaphson nr(_J, _r);
        nr.setMaxIterations(_option.max_iteration);
        nr.setAbsResidualTolerance(_option.error_tolerance);
        nr.solve(*_f_r, *_f_J, x);
    }

private:
    DISALLOW_COPY_AND_ASSIGN(NewtonRaphson);

private:
    F_R* _f_r;
    F_J* _f_J;
    MathLib::IMatrix* _J;
    MathLib::IVector* _r;
 //    T_STEP_INIT* _nl_step_init;
};



} //end

