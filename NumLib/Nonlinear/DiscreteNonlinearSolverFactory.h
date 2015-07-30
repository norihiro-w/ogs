/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 */

#pragma once

#include "NonlinearSolverOption.h"
#include "Linear.h"
#include "Picard.h"
#include "Newton.h"
#ifdef USE_PETSC
#include "PETSc/PETScSNES.h"
#endif

namespace NumLib
{

/**
 * \brief Factory class generating a discrete nonlinear solver
 */
class DiscreteNonlinearSolverFactory
{
public:
    template <class F_LINEAR, class F_R, class F_DX>
    INonlinearSolver* create(const NonlinerSolverOption &nl_option, F_LINEAR* f_l, F_R* f_r, F_DX* f_dx, MathLib::IMatrix* J, MathLib::IVector* r)
    {
        INonlinearSolver* solver = nullptr;
        switch (nl_option.solver_type)
        {
        case NonlinerSolverOption::LINEAR:
            solver = new Linear<F_LINEAR>(f_l);
            break;
        case NonlinerSolverOption::PICARD:
            solver = new Picard<F_LINEAR>(f_l);
            break;
        case NonlinerSolverOption::NEWTON:
            solver = new NewtonRaphson<F_R, F_DX>(f_r, f_dx, J, r);
            break;
#ifdef USE_PETSC
        case NonlinerSolverOption::SNES:
            solver = new PETScSNES<F_R, F_DX>(f_r, f_dx, J, r);
            break;
#endif
        default:
            ERR("*** Specified nonlinear solver not supported");
            return nullptr;
        }
        solver->setOption(nl_option);
        return solver;
    }
};

} //
