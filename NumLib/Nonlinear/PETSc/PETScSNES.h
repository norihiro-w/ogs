/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 */

#pragma once

#include "MathLib/Nonlinear/PETSc/PETScNonlinear.h"

#include "NumLib/Nonlinear/INonlinearSolver.h"

namespace NumLib
{

template <
    class F_R, 
    class F_DX
//    class T_STEP_INIT=MathLib::NRIterationStepInitializerDummy
>
class PETScSNES : public INonlinearSolver
{
public:

    PETScSNES(F_R* f_r, F_DX* f_dx, MathLib::IMatrix* J, MathLib::IVector* r)
    : _nr(J, r), _f_r(f_r), _f_dx(f_dx)
    {
    }

    virtual ~PETScSNES()
    {
        //BaseLib::releaseObject(_nl_step_init);
    };

    virtual void solve(MathLib::IVector &x_)
    {
        auto &x(static_cast<MathLib::PETScVector&>(x_));
        _nr.setMaxIterations(_option.max_iteration);
        _nr.setAbsResidualTolerance(_option.error_tolerance);
        _nr.solve(*_f_r, *_f_dx, x);
    }

private:
    DISALLOW_COPY_AND_ASSIGN(PETScSNES);

private:
    MathLib::PETScNonlinear _nr;
    F_R* _f_r;
    F_DX* _f_dx;
//    T_STEP_INIT* _nl_step_init;

};



} //end

