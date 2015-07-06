/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 */

#pragma once

#include <vector>

#include "MathLib/LinAlg/IMatrix.h"
#include "MathLib/LinAlg/ILinearSolver.h"

#include "NumLib/Function/IFunction.h"
#include "NumLib/TimeStepping/TimeStep.h"
#include "NumLib/Function/FemNodalFunction.h"

#include "FemVariable.h"
#include "FemDirichletBC.h"
#include "IFemNeumannBC.h"
#include "IElementAssembler.h"


namespace SolutionLib
{

/**
 * \brief Template class for transient linear FEM functions
 *
 * - FEM variables
 * - Local assembler
 * - Linear EQS
 * - Current time step
 * - Previous time step value
 *
 * \tparam T_ASSEMBLER
 */
class TransientFEMLinearFunction
{
public:

    /// constructor
    /// \param list_var
    /// \param assembler
    /// \param linear_eqs    Linear equation
    TransientFEMLinearFunction(
            const MeshLib::Mesh* msh,
            const std::vector<FemVariable*> &list_var,
            IElementAssembler* assembler,
            AssemblerLib::LocalToGlobalIndexMap* dofManager,
            MathLib::ILinearSolver* lsolver,
            MathLib::IMatrix* A, MathLib::IVector* x, MathLib::IVector* rhs)
    : _msh(msh), _local_assembler(assembler), _dofManager(dofManager), _linear_solver(lsolver),
      _A(A), _x(x), _rhs(rhs),
      _t_n1(0), _u_n0(0), _list_var(list_var)
    {
    }

    ///
    virtual ~TransientFEMLinearFunction() {};

    ///
    TransientFEMLinearFunction* clone() const
    {
        return new TransientFEMLinearFunction(_msh, _list_var, _local_assembler, _dofManager, _linear_solver, _A, _x, _rhs);
    }

    /// reset property
    void reset(const NumLib::TimeStep* t, const MathLib::IVector* u_n0)
    {
        this->_t_n1 = t;
        this->_u_n0 = u_n0;
    };

    /// solve linear equations discretized with FEM
    /// @param u_k0    initial guess
    /// @param u_k1    new results
    void operator()(MathLib::IVector &u);

private:
    const MeshLib::Mesh* _msh;
    IElementAssembler* _local_assembler;
    AssemblerLib::LocalToGlobalIndexMap* _dofManager;
    MathLib::ILinearSolver* _linear_solver;
    MathLib::IMatrix* _A;
	MathLib::IVector* _x;
	MathLib::IVector* _rhs;
    const NumLib::TimeStep* _t_n1;
    const MathLib::IVector* _u_n0;
    std::vector<FemVariable*> _list_var;
};


} //end
