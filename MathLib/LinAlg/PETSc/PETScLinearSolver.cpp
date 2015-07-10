/*!
   \file  PETScLinearSolver.cpp
   \brief Definition of class PETScLinearSolver, which defines a solver object
         based on PETSc routines.

   \author Wenqing Wang
   \version
   \date Nov 2011 - Mar 2014

   \copyright
   Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
               Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license
*/

#include "PETScLinearSolver.h"

#include "BaseLib/RunTime.h"
#include "PETScTools.h"

namespace MathLib
{

using boost::property_tree::ptree;

PETScLinearSolver::PETScLinearSolver(PETScMatrix &A, const boost::property_tree::ptree* option)
: _A(A), _elapsed_ctime(0.)
{
    KSPCreate(PETSC_COMM_WORLD, &_solver);
    KSPGetPC(_solver, &_pc);

    if (option)
        setOption(*option);

    const std::string prefix;
    if ( !prefix.empty() )
    {
        KSPSetOptionsPrefix(_solver, prefix.c_str());
    }
    KSPSetFromOptions(_solver);
}

void PETScLinearSolver::setOption(const boost::property_tree::ptree &option)
{
    boost::optional<ptree> ptSolver = option.get_child("LinearSolver");
     if (!ptSolver)
         return;

     boost::optional<std::string> solver_type = ptSolver->get_optional<std::string>("solver_type");
     if (solver_type) {
         _option.solver_type = _option.getSolverType(*solver_type);
     }
     boost::optional<std::string> precon_type = ptSolver->get_optional<std::string>("precon_type");
     if (precon_type) {
         _option.precon_type = _option.getPreconType(*precon_type);
     }
     boost::optional<double> error_tolerance = ptSolver->get_optional<double>("error_tolerance");
     if (error_tolerance) {
         _option.error_tolerance = *error_tolerance;
     }
     boost::optional<int> max_iteration_step = ptSolver->get_optional<int>("max_iteration_step");
     if (max_iteration_step) {
         _option.max_iterations = *max_iteration_step;
     }

}

void PETScLinearSolver::imposeKnownSolution(IMatrix &A, IVector &b, const std::vector<std::size_t> &vec_knownX_id,
        const std::vector<double> &vec_knownX_x, double /*penalty_scaling*/)
{
    applyKnownSolution(static_cast<PETScMatrix&>(A), static_cast<PETScVector&>(b), vec_knownX_id, vec_knownX_x);
}

void PETScLinearSolver::solve(MathLib::IVector &b_, MathLib::IVector &x_)
{
    Mat &A = _A.getRawMatrix();
    Vec &b = static_cast<PETScVector&>(b_).getData();
    Vec &x = static_cast<PETScVector&>(x_).getData();

    PetscPrintf(PETSC_COMM_WORLD, "------------------------------------------------------------------\n");
    PetscPrintf(PETSC_COMM_WORLD, "*** PETSc solver computation\n");

    BaseLib::RunTime wtimer;
    wtimer.start();
    	
// define TEST_MEM_PETSC
#ifdef TEST_MEM_PETSC
    PetscLogDouble mem1, mem2;
    PetscMemoryGetCurrentUsage(&mem1);
#endif

#if (PETSC_VERSION_NUMBER > 3040)
    KSPSetOperators(_solver, A, A);
#else
    KSPSetOperators(_solver, A, A, DIFFERENT_NONZERO_PATTERN);
#endif
    KSPSetType(_solver, _option.solver_type.c_str());
    PCSetType(_pc, _option.precon_type.c_str());
    KSPSetTolerances(_solver, _option.error_tolerance, PETSC_DEFAULT, PETSC_DEFAULT, _option.max_iterations);
    KSPSetFromOptions(_solver);

    KSPSolve(_solver, b, x);

//    _A.write("A.txt");
//    b_.write("b.txt");
//    x_.write("x.txt");


    //-------------------------------------------------------------------------
    // check result
    PetscPrintf(PETSC_COMM_WORLD, "solver    : %s\n", _option.solver_type.c_str());
    PetscPrintf(PETSC_COMM_WORLD, "precon    : %s\n", _option.precon_type.c_str());
    KSPConvergedReason reason;
    KSPGetConvergedReason(_solver,&reason); //CHKERRQ(ierr);

    if (reason == KSP_DIVERGED_INDEFINITE_PC)
    {
        PetscPrintf(PETSC_COMM_WORLD, "status    : Divergence because of indefinite preconditioner. Run the executable again but with -pc_factor_shift_positive_definite option.");
    }
    else if(reason == KSP_DIVERGED_ITS)
    {
        PetscPrintf(PETSC_COMM_WORLD, "status    : maximum number of iterations reached\n");
    }
    else if (reason < 0)
    {
        PetscPrintf(PETSC_COMM_WORLD, "status    : Other kind of divergence (reason=%d)\n", reason);
    }
    else
    {
        PetscPrintf(PETSC_COMM_WORLD, "status    : converged\n");
        PetscInt its;
        double res;
        KSPGetIterationNumber(_solver,&its); //CHKERRQ(ierr);
        KSPGetResidualNorm(_solver, &res);
        PetscPrintf(PETSC_COMM_WORLD, "iteration : %d/%d\n", (int)its, _option.max_iterations);
        PetscPrintf(PETSC_COMM_WORLD, "residual  : %e\n", res);
    }
    PetscPrintf(PETSC_COMM_WORLD, "------------------------------------------------------------------\n");

#ifdef TEST_MEM_PETSC
    PetscMemoryGetCurrentUsage(&mem2);
    PetscPrintf(PETSC_COMM_WORLD, "###Memory usage by solver. Before :%f After:%f Increase:%d\n", mem1, mem2, (int)(mem2 - mem1));
#endif

    _elapsed_ctime += wtimer.elapsed();

    //return converged;
}

} //end of namespace

