/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 */

#pragma once

#include <vector>

#include <logog/include/logog.hpp>


#include "MathLib/LinAlg/LinAlgBuilder.h"

#include "MeshLib/Mesh.h"
#include "MeshLib/MeshSubsets.h"

#include "AssemblerLib/LocalToGlobalIndexMap.h"
#include "AssemblerLib/MapTools.h"

#include "NumLib/Function/FemNodalFunction.h"
#include "NumLib/Nonlinear/TemplateDiscreteNonlinearSolver.h"
#include "NumLib/Nonlinear/DiscreteNonlinearSolverFactory.h"

#include "SolutionLib/Core/AbstractTimeSteppingAlgorithm.h"
#include "FemDirichletBC.h"
#include "FemNeumannBC.h"
#include "FemIVBVProblem.h"


namespace SolutionLib
{
class TransientFEMLinearFunction;
class TransientFEMResidualFunction;
class TransientFEMJacobianFunction;

/**
 * \brief Solution algorithm for linear transient problems using FEM with single time stepping method
 *
 * - previous and current time step values
 * - ST values
 * - Linear equation
 * - Residual equation
 * - Dx equation
 *
 */
class SingleStepFEM
    : public AbstractTimeSteppingAlgorithm
{
public:
    typedef NumLib::TemplateDiscreteNonlinearSolver<TransientFEMLinearFunction, TransientFEMResidualFunction, TransientFEMJacobianFunction, NumLib::DiscreteNonlinearSolverFactory> NonlinearSolverType;
    typedef typename NumLib::FemNodalFunctionScalar MyNodalFunctionScalar;

    /// constructor
    /// - initialize solution vectors
    /// - set up DoFs and equation index
    /// - prepare linear equation and solver
    /// - prepare linear functions
    SingleStepFEM(MeshLib::Mesh* msh,
            FemIVBVProblem* problem,
            MathLib::LinAlgLibType lsLibType);

    virtual ~SingleStepFEM();

    /// solve 
    int solveTimeStep(const NumLib::TimeStep &t_n1);

    /// get the current solution
    MyNodalFunctionScalar* getCurrentSolution(size_t var_id) { return _vec_u_n1[var_id]; }

    ///
    FemIVBVProblem* getProblem() {return _problem;};

    /// get a linear solver
    MathLib::ILinearSolver* getLinearEquationSolver() { return _linear_solver; };

    /// get a nonlinear solver
    NonlinearSolverType* getNonlinearSolver() { return _f_nonlinear;};

    AssemblerLib::LocalToGlobalIndexMap* getDofEquationIdTable() {return _dofManager;};

    ///
    virtual void accept(const NumLib::TimeStep &t)
    {
        AbstractTimeSteppingAlgorithm::accept(t);
        *_x_n0 = *_x_n1; //copy current value to previous value
    };

private:
    DISALLOW_COPY_AND_ASSIGN(SingleStepFEM);


private:
    MeshLib::Mesh* _msh;
    AssemblerLib::LocalToGlobalIndexMap* _dofManager;
    FemIVBVProblem* _problem;
    MathLib::ILinearSolver* _linear_solver;
    MathLib::IMatrix* _A;
    MathLib::IVector* _b;
    MathLib::IVector* _x;
    std::vector<MyNodalFunctionScalar*> _vec_u_n1;
    TransientFEMLinearFunction* _f_linear;
    TransientFEMResidualFunction* _f_r;
    TransientFEMJacobianFunction* _f_J;
    NonlinearSolverType* _f_nonlinear;
    MathLib::IVector *_x_n0;
    MathLib::IVector *_x_n1_0;
    MathLib::IVector *_x_n1;
    MathLib::IVector *_x_st;
};

}
