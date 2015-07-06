/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 */

#pragma once

#include <vector>

namespace MathLib
{
class IVector;
class IMatrix;
class ILinearSolver;
}

namespace MeshLib
{
class Mesh;
}

namespace AssemblerLib
{
class LocalToGlobalIndexMap;
}

namespace NumLib
{
class TimeStep;
}

namespace SolutionLib
{
class FemVariable;
class IElementAssembler;

/**
 *  Form a Jacobian matrix
 */
class TransientFEMJacobianFunction
{
public:
    TransientFEMJacobianFunction(
        const MeshLib::Mesh* msh,
        const std::vector<FemVariable*> &list_var,
        IElementAssembler* asssembler,
        AssemblerLib::LocalToGlobalIndexMap* dofManager)
        : _msh(msh), _local_assembler(asssembler), _dofManager(dofManager),
          _t_n1(nullptr), _u_n0(nullptr), _list_var(list_var)
    {
    }

    virtual ~TransientFEMJacobianFunction() {}

    TransientFEMJacobianFunction* clone() const
    {
        return new TransientFEMJacobianFunction(_msh, _list_var, _local_assembler, _dofManager);
    }

    /// reset property
    void reset(const NumLib::TimeStep* t, const MathLib::IVector* u_n0)
    {
        this->_t_n1 = t;
        this->_u_n0 = u_n0;
    }

    /// form Jocobian
    /// @param u_n1    current solution
    /// @param J       Jacobian matrix to be formed
    void operator()(const MathLib::IVector &u_n1, MathLib::IMatrix &J);

private:
    const MeshLib::Mesh* _msh;
    IElementAssembler *_local_assembler;
    AssemblerLib::LocalToGlobalIndexMap* _dofManager;
    const NumLib::TimeStep* _t_n1;
    const MathLib::IVector* _u_n0;
    std::vector<FemVariable*> _list_var;
};


} //end
