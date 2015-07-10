/*!
   \file  PETScLinearSolver.h
   \brief Declaration of class PETScLinearSolver, which defines a solver object
         based on PETSc routines.

   \author Wenqing Wang
   \version
   \date Nov 2011 - Sep 2013

   \copyright
   Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
               Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license
*/

#ifndef PETSCLINEARSOLVER_H_
#define PETSCLINEARSOLVER_H_

#include <string>

#include <logog/include/logog.hpp>
#include <petscksp.h>

#include "MathLib/LinAlg/ILinearSolver.h"

#include "PETScMatrix.h"
#include "PETScVector.h"
#include "PETScOption.h"

namespace MathLib
{

/*!
     A class of linear solver based on PETSc rountines.

     All options for KSP and PC must given in the command line.
*/
class PETScLinearSolver : public ILinearSolver
{
    public:

        PETScLinearSolver(PETScMatrix &A, const boost::property_tree::ptree* option=nullptr);

        ~PETScLinearSolver()
        {
            INFO("Time elapsed in PETSc ksp solver for equation: %g s.\n",
                 _elapsed_ctime);

            KSPDestroy(&_solver);
        }

        LinAlgLibType getLinAlgLibType() const { return LinAlgLibType::PETSc; }

        /**
         * configure linear solvers
         * @param option
         */
        void setOption(const boost::property_tree::ptree &option);

        /*!
            Solve a system of equations.
            \param b The right hand side of the equations.
            \param x The solutions to be solved.
            \return  true: converged, false: diverged.
        */
        void solve(MathLib::IVector &b, MathLib::IVector &x);

        /// apply prescribed values to a system of linear equations
        void imposeKnownSolution(IMatrix &A, IVector &b, const std::vector<std::size_t> &vec_knownX_id,
                const std::vector<double> &vec_knownX_x, double penalty_scaling = 1e+10);

        /// Get number of iterations.
        PetscInt getNIterations() const
        {
            PetscInt its = 0;
            KSPGetIterationNumber(_solver, &its);
            return its;
        }

        /// Get elapsed wall clock time.
        double getElapsedTime() const
        {
            return _elapsed_ctime;
        }

    private:
        /// Matrix, kept as a member only for solving successive linear equation
        /// that preconditioner matrix may vary.
        PETScMatrix &_A;

        KSP _solver; ///< Solver type.
        PC _pc;      ///< Preconditioner type.

        PETScOption _option;
        double _elapsed_ctime; ///< Clock time
};

} // end namespace
#endif

