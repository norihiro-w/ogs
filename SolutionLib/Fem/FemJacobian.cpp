/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 */

#include "FemJacobian.h"

#include "MathLib/LinAlg/ILinearSolver.h"
#include "MathLib/LinAlg/IMatrix.h"
#include "MathLib/LinAlg/IVector.h"

#include "MeshLib/Elements/Element.h"

#include "AssemblerLib/LocalToGlobalIndexMap.h"
#include "AssemblerLib/MeshComponentMap.h"
#include "AssemblerLib/MapTools.h"

#include "NumLib/Function/IFunction.h"
#include "NumLib/TimeStepping/TimeStep.h"
#include "NumLib/Function/FemNodalFunction.h"

#include "IElementAssembler.h"
#include "FemDirichletBC.h"
#include "FemNeumannBC.h"
#include "FemVariable.h"


namespace SolutionLib
{

void TransientFEMJacobianFunction::operator()(const MathLib::IVector &u_n1, MathLib::IMatrix &J)
{
    // input, output
    const NumLib::TimeStep &t_n1 = *this->_t_n1;
    const MathLib::IVector &u_n = *this->_u_n0;
    const size_t msh_id = _msh->getID();

    J.setZero();

    // assembly
    for (auto e : _msh->getElements()) {
        // get dof map
        auto rowColIndeces = (*_dofManager)[e->getID()];

        // local assembly
        MathLib::LocalMatrix localA = MathLib::LocalMatrix::Zero(rowColIndeces.rows.size(), rowColIndeces.rows.size());
        // previous and current results
        auto local_u_n1 = AssemblerLib::getLocalVector(rowColIndeces.rows, u_n1);
        auto local_u_n = AssemblerLib::getLocalVector(rowColIndeces.columns, u_n);

        // create local DoF table
        std::vector<MeshLib::Node*> vec_items;
        for (std::size_t i=0; i<e->getNNodes(); i++)
            vec_items.push_back(const_cast<MeshLib::Node*>(e->getNode(i)));
        MeshLib::MeshSubset subset(*_msh, &vec_items);
        MeshLib::MeshSubsets ss(&subset);
        std::vector<MeshLib::MeshSubsets*> mesh_subsets(1, &ss);
        AssemblerLib::LocalToGlobalIndexMap localDofMap(mesh_subsets);

        _local_assembler->reset(*e, localDofMap);
        _local_assembler->jacobian(t_n1, local_u_n1, local_u_n, localA);

        // update global
        J.add(rowColIndeces.rows, rowColIndeces.columns, localA);
    }
    J.assemble();

    // set residual to zero for Dirichlet BC
    std::vector<size_t> list_bc1_eqs_id;
    std::vector<double> list_bc1_val;
    for (size_t i=0; i<_list_var.size(); i++) {
        FemVariable* var = _list_var[i];
        std::vector<size_t> var_bc_id;
        std::vector<double> var_bc_val;
        for (size_t j=0; j<var->getNumberOfDirichletBC(); j++) {
            FemDirichletBC* bc1 = var->getDirichletBC(j);
            bc1->setup(var->getCurrentOrder());
            AssemblerLib::toGlobalIDValues(*_dofManager, i, msh_id, bc1->getListOfBCNodes(), bc1->getListOfBCValues(), list_bc1_eqs_id, list_bc1_val);
        }
    }
    J.zeroRows(list_bc1_eqs_id);

    J.assemble();
}

} //end
