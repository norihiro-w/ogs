/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 */

#include "FemLinearEQS.h"

#include "AssemblerLib/MapTools.h"

namespace SolutionLib
{

void TransientFEMLinearFunction::operator()(MathLib::IVector &u_k1)
{
    // input, output
    const NumLib::TimeStep &t_n1 = *this->_t_n1;
    const MathLib::IVector &u_n = *this->_u_n0;

    _A->setZero();
    (*_rhs) = 0.0;

    // setup BC
    for (size_t i=0; i<_list_var.size(); i++) {
    	FemVariable* var = _list_var[i];
        for (size_t j=0; j<var->getNumberOfNeumannBC(); j++) {
            IFemNeumannBC* st = var->getNeumannBC(j);
            st->initCurrentTime(t_n1.current());
            st->setup(var->getCurrentOrder());
        }
    }

    // assembly
    MathLib::LocalMatrix localA;
    MathLib::LocalVector localRhs;
    MathLib::LocalVector local_u_n1, local_u_n;
    for (auto e : _msh->getElements())
    {
        // get dof map
        const auto rowColIndeces((*_dofManager)[e->getID()]);

        // previous and current results
        local_u_n1.resize(rowColIndeces.columns.size());
        local_u_n.resize(rowColIndeces.columns.size());
        AssemblerLib::getLocalVector(rowColIndeces.columns, u_k1, local_u_n1);
        AssemblerLib::getLocalVector(rowColIndeces.columns, u_n, local_u_n);

#if 0
        // create a local DoF table
        std::vector<MeshLib::Node*> vec_items;
        for (std::size_t i=0; i<e->getNNodes(); i++)
            vec_items.push_back(const_cast<MeshLib::Node*>(e->getNode(i)));
        MeshLib::MeshSubset subset(*_msh, &vec_items);
        MeshLib::MeshSubsets ss(&subset);
        std::vector<MeshLib::MeshSubsets*> mesh_subsets(1, &ss);
        AssemblerLib::LocalToGlobalIndexMap localDofMap(mesh_subsets, AssemblerLib::ComponentOrder::BY_COMPONENT, false);
#endif
        // local assembly
        localA.resize(rowColIndeces.rows.size(), rowColIndeces.rows.size());
        localRhs.resize(rowColIndeces.rows.size());
        localA.setZero();
        localRhs.setZero();
//        MathLib::LocalMatrix localA(rowColIndeces.rows.size(), rowColIndeces.rows.size());
//        MathLib::LocalVector localRhs(rowColIndeces.rows.size());
        _local_assembler->reset(*e);
        _local_assembler->linear(t_n1, local_u_n1, local_u_n, localA, localRhs);

        // update global
        _A->add(rowColIndeces.rows, rowColIndeces.columns, localA);
        _rhs->add(rowColIndeces.rows, localRhs);

//        std::stringstream sst; sst << rowColIndeces;
//        sst << "A=\n" <<  localA;
//        sst << "\nb=\n" << localRhs;
//        INFO("e %d:\%s", e->getID(), sst.str().data());
    }
    _A->assemble();
    _rhs->assemble();

    //apply BC1,2
    for (size_t i=0; i<_list_var.size(); i++) {
    	FemVariable* var = _list_var[i];
        std::vector<size_t> var_bc_id;
        std::vector<double> var_bc_val;
        for (size_t j=0; j<var->getNumberOfDirichletBC(); j++) {
            FemDirichletBC* bc1 = var->getDirichletBC(j);
            bc1->setup(var->getCurrentOrder());
            AssemblerLib::toGlobalIDValues(*_dofManager, i, _msh->getID(), bc1->getListOfBCNodes(), bc1->getListOfBCValues(), var_bc_id, var_bc_val);
        }
        _linear_solver->imposeKnownSolution(*_A, *_rhs, var_bc_id, var_bc_val);

        for (size_t j=0; j<var->getNumberOfNeumannBC(); j++) {
            IFemNeumannBC* bc2 = var->getNeumannBC(j);
            std::vector<size_t> var_bc2_id;
            std::vector<double> var_bc2_val;
            AssemblerLib::toGlobalIDValues(*_dofManager, i, _msh->getID(), bc2->getListOfBCNodes(), bc2->getListOfBCValues(), var_bc2_id, var_bc2_val);
            _rhs->add(var_bc2_id, var_bc2_val, -1.0);
        }
    }
    _A->assemble();
    _rhs->assemble();

//#define OUT_LEQS
#ifdef OUT_LEQS
    _A->assemble();
    _A->write("A.txt");
    _rhs->write("b.txt");
#endif

    // solve
    *_x = u_k1;
    _linear_solver->solve(*_rhs, *_x);
    u_k1 = *_x;

#ifdef OUT_LEQS
    _x->write("x.txt");
#endif
}

} //end
