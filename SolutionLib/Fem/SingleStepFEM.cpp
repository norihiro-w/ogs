/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 */

#include "SingleStepFEM.h"

#include "BaseLib/DebugTools.h"

#include "MathLib/LinAlg/ILinearSolver.h"
#include "MathLib/LinAlg/IMatrix.h"
#include "MathLib/LinAlg/IVector.h"

#include "AssemblerLib/LocalToGlobalIndexMap.h"
#include "AssemblerLib/VecMatBuilder.h"

#include "FemLinearEQS.h"
#include "FemResidualEQS.h"
#include "FemJacobian.h"

namespace SolutionLib
{

SingleStepFEM::SingleStepFEM(
		MeshLib::Mesh* msh, FemIVBVProblem* problem, MathLib::LinAlgLibType lsLibType)
    : AbstractTimeSteppingAlgorithm(*problem->getTimeSteppingFunction()),
      _msh(msh), _problem(problem)
{
    INFO("->setting up a solution algorithm SingleStepFEM");

    const size_t n_var = problem->getNumberOfVariables();

    // create dof map
    INFOa("->constructing a DoF table");
    std::vector<MeshLib::MeshSubsets*> components;
    for (size_t i=0; i<n_var; i++) {
        FemVariable* var = problem->getVariable(i);
        size_t n_dof_per_var = msh->getNNodes(); //(var->getCurrentOrder());
        components.push_back(new MeshLib::MeshSubsets(new MeshLib::MeshSubset(*msh, &msh->getNodes())));
        INFOa("* Variable %d: name=%s, order=%d, n_dof=%d", i, var->getName().c_str(), var->getCurrentOrder(), n_dof_per_var);
    }
    _dofManager = new AssemblerLib::LocalToGlobalIndexMap(components, AssemblerLib::ComponentOrder::BY_COMPONENT);
    const size_t n_total_dofs = _dofManager->dofSize();
    INFOa("* Total number of DoFs = %d", n_total_dofs);
//    std::stringstream ss;
//    ss << _dofManager->getMeshComponentMap();
//    INFO("%s", ss.str().c_str())

    // setup IC
    _vec_u_n1.resize(n_var, nullptr);
    for (size_t i=0; i<n_var; i++) {
    	FemVariable* femVar = problem->getVariable(i);
        FemIC* femIC = femVar->getIC();
        MyNodalFunctionScalar* u0 = new MyNodalFunctionScalar(lsLibType); //TODO ddc, different from solution vectors
        u0->initialize(_msh, femVar->getCurrentOrder());
        assert (femVar->getFeObjectContainer() != NULL);
        u0->setFeObjectContainer(femVar->getFeObjectContainer());
        femIC->setup(*u0->getNodalValues());
        _vec_u_n1[i] = u0; //u0->clone()
    }

    // initialize vectors for solution, ST
    _x_n0 = AssemblerLib::VecMatBuilder::generateVector(lsLibType, _dofManager->getMeshComponentMap());
    _x_n1 = _x_n0->duplicate();
    _x_n1_0 = _x_n0->duplicate();
    _x_st = _x_n0->duplicate();

    // copy values of each variable to one solution vector
    for (size_t i=0; i<n_var; i++) {
        MathLib::IVector* vec_var = _vec_u_n1[i]->getNodalValues();
        AssemblerLib::setGlobalVector(*_dofManager, i, msh->getID(), *vec_var, *_x_n0);
    }

    // create linear equation systems
    _x = _x_n0->duplicate();
    _b = _x->duplicate();
    _A = AssemblerLib::VecMatBuilder::generateMatrix(lsLibType, _dofManager->getMeshComponentMap());
    _linear_solver = MathLib::LinAlgBuilder::generateLinearSolver(lsLibType, _A);

    // setup functions
    std::vector<FemVariable*> list_var(n_var);
    for (size_t i=0; i<n_var; i++) list_var[i] = problem->getVariable(i);
    _f_linear = new TransientFEMLinearFunction(msh, list_var, problem->getAssembler(), _dofManager, _linear_solver, _A, _x, _b);
    _f_r = new TransientFEMResidualFunction(msh, list_var, _dofManager, problem->getAssembler());
    _f_J = new TransientFEMJacobianFunction(msh, list_var, problem->getAssembler(), _dofManager);

    // setup nonlinear solver
    _f_nonlinear = new NonlinearSolverType(_f_linear, _f_r, _f_J, _A, _b);
};

SingleStepFEM::~SingleStepFEM()
{
    delete _linear_solver;
    delete _x_n0;
    delete _x_n1;
    delete _x_n1_0;
    delete _x_st;
    delete _x;
    delete _b;
    delete _A;
    delete _f_linear;
    delete _f_r;
    delete _f_J;
    delete _f_nonlinear;
    delete _dofManager;
}

int SingleStepFEM::solveTimeStep(const NumLib::TimeStep &t_n1)
{
    BaseLib::MPIEnvironment mpi;
    // time step
//    double dt = t_n1.dt(); // - AbstractTimeSteppingAlgorithm::getTimeStepFunction()->getPrevious();
    NumLib::TimeStep this_t_n1(t_n1);

    const size_t n_var = _problem->getNumberOfVariables();
    const size_t msh_id = _msh->getID();
    // bc1
    std::vector<size_t> list_bc1_eqs_id;
    std::vector<double> list_bc1_val;
    for (size_t i_var=0; i_var<n_var; i_var++) {
        FemVariable* var = _problem->getVariable(i_var);
        for (size_t i=0; i<var->getNumberOfDirichletBC(); i++) {
            SolutionLib::FemDirichletBC *bc1 = var->getDirichletBC(i);
            bc1->setup(var->getCurrentOrder());
            INFO("-> setting 1st BC: %s, %d nodes", var->getName().data(), bc1->getListOfBCNodes().size());
            //std::stringstream ss; bc1->write(ss); INFO("%s", ss.str().data());
            std::vector<size_t> &list_bc_nodes = bc1->getListOfBCNodes();
            std::vector<double> &list_bc_values = bc1->getListOfBCValues();
            AssemblerLib::toGlobalIDValues(*_dofManager, i_var, msh_id, list_bc_nodes, list_bc_values, list_bc1_eqs_id, list_bc1_val);
        }
    }

    // st
    std::vector<size_t> list_st_eqs_id;
    std::vector<double> list_st_val;
    for (size_t i_var=0; i_var<n_var; i_var++) {
    	FemVariable* var = _problem->getVariable(i_var);
        for (size_t i=0; i<var->getNumberOfNeumannBC(); i++) {
            SolutionLib::IFemNeumannBC *bc2 = var->getNeumannBC(i);
            bc2->initCurrentTime(this_t_n1.current());
            bc2->setup(var->getCurrentOrder());
            INFO("-> setting 2nd BC: %s, %d nodes", var->getName().data(), bc2->getListOfBCNodes().size());
            std::vector<size_t> &list_bc_nodes = bc2->getListOfBCNodes();
            std::vector<double> &list_bc_values = bc2->getListOfBCValues();
//            std::cout << list_bc_nodes << std::endl;
//            std::cout << list_bc_values << std::endl;
            AssemblerLib::toGlobalIDValues(*_dofManager, i_var, msh_id, list_bc_nodes, list_bc_values, list_st_eqs_id, list_st_val);
        }
    }
    (*_x_st) = .0;
    for (size_t i=0; i<list_st_eqs_id.size(); i++) {
        _x_st->set(list_st_eqs_id[i], list_st_val[i]);
    }
    _x_st->assemble();

    // setup functions
    _f_linear->reset(&t_n1, _x_n0);
    _f_r->reset(&t_n1, _x_n0, _x_st);
    _f_J->reset(&t_n1, _x_n0);

    // initial guess
    if (mpi.root()) INFO("->use previous time step value as initial guess");
    *_x_n1_0 = *_x_n0;
    for (size_t i=0; i<list_bc1_eqs_id.size(); i++) {
    	_x_n1_0->set(list_bc1_eqs_id[i], list_bc1_val[i]);
    }

    // solve
    *_x_n1 = *_x_n1_0;
    _f_nonlinear->solve(*_x_n1);
    _x_n1->assemble();

    // distribute solution vector to local vector for each variable
    for (size_t i=0; i<n_var; i++) {
        AssemblerLib::setLocalVector(*_dofManager, i, msh_id, *_x_n1, *_vec_u_n1[i]->getNodalValues());
        _vec_u_n1[i]->getNodalValues()->assemble();
    }

    return 0;
}

}
