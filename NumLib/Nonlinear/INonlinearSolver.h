/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file INonlinearSolver.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include "MathLib/DataType.h"

#include "NonlinearSolverOption.h"

namespace NumLib
{

/**
 * \brief Interface to nonlinear solvers
 */
class INonlinearSolver
{
public:
    virtual ~INonlinearSolver() {};

    virtual void solve(MathLib::IVector &x_new) = 0;

    NonlinerSolverOption& getOption() {return _option;};
    void setOption(const NonlinerSolverOption &option) {_option = option;};
protected:
    NonlinerSolverOption _option;
};

} //end

