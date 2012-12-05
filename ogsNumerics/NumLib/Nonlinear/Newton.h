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
#include "MathLib/Nonlinear/NRIterationStepInitializerDummy.h"
#include "DiscreteLib/Core/IDiscreteSystem.h"

#include "INonlinearSolver.h"

namespace NumLib
{

/**
 * \brief Newton-Raphson solver for discrete systems
 */
template <
    class T_DIS_SYS, 
    class F_R, 
    class F_DX, 
    class T_STEP_INIT=MathLib::NRIterationStepInitializerDummy 
>
class NewtonRaphson : public INonlinearSolver
{
public:
    typedef T_DIS_SYS MyDiscreteSystem;

    NewtonRaphson(MyDiscreteSystem* dis_sys, F_R* f_r, F_DX* f_dx)
    : _dis_sys(dis_sys), _f_r(f_r), _f_dx(f_dx), _nl_step_init(NULL)
    {
        _r = _dx = 0;
    };

    NewtonRaphson(MyDiscreteSystem* dis_sys, F_R* f_r, F_DX* f_dx, T_STEP_INIT* nl_step_init)
    : _dis_sys(dis_sys), _f_r(f_r), _f_dx(f_dx), _nl_step_init(nl_step_init)
    {
        _r = _dx = 0;
    };

    virtual ~NewtonRaphson()
    {
        _dis_sys->deleteVector(_r);
        _dis_sys->deleteVector(_dx);
        BaseLib::releaseObject(_nl_step_init);
    };

    virtual void solve(const VectorType &x_0, VectorType &x_new)
    {
        if (_r==0) {
            _r = _dis_sys->template createVector<double>(x_0.size());
            _dx = _dis_sys->template createVector<double>(x_0.size());
        }
        MathLib::NRCheckConvergence<VectorType, MathLib::NRErrorAbsResMNormOrRelDxMNorm> check(_option.error_tolerance); 
        // MathLib::NRCheckConvergence<VectorType, MathLib::NRErrorAbsResMNorm> check(_option.error_tolerance);  // HS
        MathLib::NewtonRaphsonMethod nr;
        nr.solve(*_f_r, *_f_dx, x_0, x_new, *_r, *_dx, _option.max_iteration, &check, _nl_step_init);
    }

private:
    DISALLOW_COPY_AND_ASSIGN(NewtonRaphson);

private:
    MyDiscreteSystem* _dis_sys;
    F_R* _f_r;
    F_DX* _f_dx;
    VectorType* _r;
    VectorType* _dx;
    T_STEP_INIT* _nl_step_init;
};



} //end

