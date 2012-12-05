/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file SingleStepKinReduction.h
 *
 * Created on 2012-09-20 by Haibing Shao
 */

#ifndef SINGLE_STEP_KIN_REDUCTION_H
#define SINGLE_STEP_KIN_REDUCTION_H

#include <vector>
#include <map>
#include "logog.hpp"

#include "BaseLib/CodingTools.h"
#include "MeshLib/Core/IMesh.h"
#include "DiscreteLib/Core/IDiscreteSystem.h"
#include "DiscreteLib/Core/IDiscreteLinearEquation.h"
#include "DiscreteLib/Utils/SparsityBuilderFromNodeConnectivity.h"
#include "FemLib/Function/FemFunction.h"
#include "NumLib/Nonlinear/TemplateDiscreteNonlinearSolver.h"

#include "SolutionLib/DataType.h"
#include "SolutionLib/Core/AbstractTimeSteppingAlgorithm.h"
#include "SolutionLib/Fem/TemplateFemEquation.h"
#include "SolutionLib/Fem/FemDirichletBC.h"
#include "SolutionLib/Fem/FemNeumannBC.h"

#include "ReductionKinNodeInfo.h"

class ReductionKinNodeInfo; 

namespace SolutionLib
{

/**
 * \brief Solution algorithm for kinetic reaction reduction problems using FEM with single time stepping method
 *
 * - previous and current time step values
 * - ST values
 * - Linear equation
 * - Residual equation
 * - Dx equation
 * - Linear EQS
 * - Linear solver
 *
 * \tparam T_USER_FEM_PROBLEM      the FEM problem class
 * \tparam T_LINEAR_SOLVER         Linear equation solver class
 */
template <
	class T_USER_FUNCTION_DATA, 
    class T_USER_FEM_PROBLEM,
    class T_USER_LINEAR_PROBLEM,
    class T_USER_LINEAR_SOLUTION,
    class T_USER_NON_LINEAR_PROBLEM, 
	class T_USER_NON_LINEAR_SOLUTION 
    >
class SingleStepKinReduction
    : public AbstractTimeSteppingAlgorithm
{
public:

	typedef T_USER_FUNCTION_DATA       UserFunctionData; 
    typedef T_USER_FEM_PROBLEM         UserFemProblem;
    typedef T_USER_LINEAR_PROBLEM      UserLinearProblem; 
	typedef T_USER_LINEAR_SOLUTION     UserLinearSolution;
    typedef T_USER_NON_LINEAR_PROBLEM  UserNonLinearProblem; 
	typedef T_USER_NON_LINEAR_SOLUTION UserNonLinearSolution; 

	typedef typename UserFemProblem::MyDiscreteSystem MyDiscreteSystem;
    typedef typename UserFemProblem::MyVariable MyVariable;
    typedef typename FemLib::FemNodalFunctionScalar<MyDiscreteSystem>::type MyNodalFunctionScalar;

    /** 
      * constructor, including the following tasks
      * - initialize solution vectors
      * - set up DoFs and equation index
      * - prepare linear and nonlinear equations
      * - prepare linear and nonlinear solutions
      */
    SingleStepKinReduction(MyDiscreteSystem*                dis, 
		                   UserFemProblem*                  problem, 
						   UserFunctionData*                function_data, 
                           std::vector<UserLinearProblem*>  linear_problem,
						   std::vector<UserLinearSolution*> linear_solutions,
                           UserNonLinearProblem*            non_linear_problem,
						   UserNonLinearSolution*           non_linear_solution);

    /**
      * destructor, reclaiming the memory
      */ 
    virtual ~SingleStepKinReduction()
    {
        // size_t i; 
        _discrete_system->deleteVector(_x_n0);
        _discrete_system->deleteVector(_x_n1);
        _discrete_system->deleteVector(_x_n1_0);
        _discrete_system->deleteVector(_x_st);

        BaseLib::releaseObject(_feObjects);
        BaseLib::releaseObject(_non_linear_problem);
        BaseLib::releaseObject(_nlin_solution);
        BaseLib::releaseObject(_function_data);

        BaseLib::releaseObjectsInStdVector( _linear_problem );
        BaseLib::releaseObjectsInStdVector( _lin_solutions );
        BaseLib::releaseObjectsInStdVector( _vec_u_n1 ); 

        BaseLib::releaseObjectsInStdMap( _bc_info ); 
    }

    /**
      * solve the time step
      */ 
    int solveTimeStep(const NumLib::TimeStep &t_n1);
    
    /**
	  * called when this solution is accepted
	  */
    virtual void accept(const NumLib::TimeStep &t)
    {
        // this solution itself
        AbstractTimeSteppingAlgorithm::accept(t);
        // call all linear solutions
        for (size_t i=0; i < _lin_solutions.size(); i++ )
            _lin_solutions[i]->accept(t); 
        // the non-linear solution
        _nlin_solution->accept(t); 
    }

    /** 
      * get the corresponding solution pointer
      */ 
    MyNodalFunctionScalar* getCurrentSolution(size_t var_id) { return _vec_u_n1[var_id]; }

    /**
      * get the pointer to the problem class
      */ 
    UserFemProblem* getProblem() {return _problem;};

    /** 
      * get the domain of freedom equation ID table
      */ 
    DiscreteLib::DofEquationIdTable* getDofEquationIdTable() {return &_dofManager;};

    /**
      * tell whether a particular node is a boundary node
      */ 
	bool isBCNode(size_t node_idx); 

private:
    DISALLOW_COPY_AND_ASSIGN(SingleStepKinReduction);

private:
    /**
      * DOF table
      */ 
    DiscreteLib::DofEquationIdTable _dofManager;

    /**
      * pointer to the problem
      */ 
    UserFemProblem* _problem;

    /**
      * pointer to discretization
      */ 
    MyDiscreteSystem *_discrete_system;


    std::vector<MyNodalFunctionScalar*> _vec_u_n1;
    SolutionVector *_x_n0;
    SolutionVector *_x_n1_0;
    SolutionVector *_x_n1;
    SolutionVector *_x_st;

    /**
      * pointer to FEM object
      */ 
    FemLib::LagrangianFeObjectContainer* _feObjects;
    
    /**
      * pointer to linear problems
      */ 
    std::vector<UserLinearProblem*>  _linear_problem;

    /**
      * pointer to linear solutions
      */ 
	std::vector<UserLinearSolution*> _lin_solutions; 

    /**
      * pointer to non-linear problem
      */ 
	UserNonLinearProblem*            _non_linear_problem;

    /**
      * pointer to non-linear solution
      */ 
    UserNonLinearSolution*           _nlin_solution; 

    /**
      * pointer to the function concentrations
      */ 
	UserFunctionData*                _function_data; 

    /**
      * boundary condition values
      */ 
	std::map<size_t, ReductionKinNodeInfo*> _bc_info; 
};

template <
	class T_USER_FUNCTION_DATA, 
    class T_USER_FEM_PROBLEM,
    class T_USER_LINEAR_PROBLEM,
    class T_USER_LINEAR_SOLUTION,
    class T_USER_NON_LINEAR_PROBLEM, 
	class T_USER_NON_LINEAR_SOLUTION 
    >
SingleStepKinReduction<T_USER_FUNCTION_DATA, T_USER_FEM_PROBLEM, T_USER_LINEAR_PROBLEM, T_USER_LINEAR_SOLUTION, T_USER_NON_LINEAR_PROBLEM, T_USER_NON_LINEAR_SOLUTION>
    ::SingleStepKinReduction(MyDiscreteSystem*                dis, 
                             UserFemProblem*                  problem, 
							 UserFunctionData*                function_data, 
                             std::vector<UserLinearProblem*>  linear_problem,
						     std::vector<UserLinearSolution*> linear_solutions,
                             UserNonLinearProblem*            non_linear_problem,
						     UserNonLinearSolution*           non_linear_solution)
    : AbstractTimeSteppingAlgorithm(*problem->getTimeSteppingFunction()),
      _problem(problem), _discrete_system(dis), _function_data(function_data), 
      _linear_problem(linear_problem), _lin_solutions(linear_solutions), 
      _non_linear_problem(non_linear_problem), _nlin_solution(non_linear_solution)
{
    INFO("->Setting up a solution algorithm SingleStepKinReduction");

    size_t i, j, i_var; 
    const size_t n_var = problem->getNumberOfVariables();
    MeshLib::IMesh* msh = dis->getMesh();

    // create dof map
    for (size_t i=0; i<n_var; i++) {
        MyVariable* var = problem->getVariable(i);
        size_t n_dof_per_var = msh->getNumberOfNodes(var->getCurrentOrder());
        _dofManager.addVariableDoFs(msh->getID(), 0, n_dof_per_var);
        INFO("* Variable %d: name=%s, order=%d, n_dof=%d", i, var->getName().c_str(), var->getCurrentOrder(), n_dof_per_var);
    }
    _dofManager.construct();
    const size_t n_total_dofs = _dofManager.getTotalNumberOfActiveDoFs();
    INFO("* Total number of DoFs = %d", n_total_dofs);

    _feObjects = new FemLib::LagrangianFeObjectContainer(*msh);

    // set up initial condition
    _vec_u_n1.resize(n_var, 0);
    for (size_t i=0; i<n_var; i++) {
        MyVariable* femVar = problem->getVariable(i);
        FemIC* femIC = femVar->getIC();
        MyNodalFunctionScalar* u0 = new MyNodalFunctionScalar();
        u0->initialize(*dis, femVar->getCurrentOrder(), 0);
        //u0->setFeObjectContainer(_feObjects);
        femIC->setup(*u0);
        u0->setFeObjectContainer(_feObjects);
        _vec_u_n1[i] = u0; //u0->clone()
    }

    // initialize vectors for solution and ST
    _x_n0 = dis->template createVector<double>(n_total_dofs);
    _x_n1 = dis->template createVector<double>(n_total_dofs);
    _x_n1_0 = dis->template createVector<double>(n_total_dofs);
    _x_st = dis->template createVector<double>(n_total_dofs);

    // copy values of each variable to one solution vector
    for ( i=0; i<n_var; i++) {
        SolutionVector* vec_var = _vec_u_n1[i]->getDiscreteData();
        DiscreteLib::setGlobalVector(_dofManager, i, msh->getID(), *vec_var, *_x_n0);
    }

    // setup functions
    std::vector<MyVariable*> list_var(n_var);
    for ( i=0; i<n_var; i++) 
		list_var[i] = problem->getVariable(i);

	// getting the boundary conditions of concentrations for all components, 
    const size_t msh_id = _discrete_system->getMesh()->getID();
    std::vector<size_t> list_bc1_eqs_id;
    std::vector<double> list_bc1_val;
    for ( i_var=0; i_var < n_var; i_var++ ) {
        MyVariable* var = _problem->getVariable(i_var);
        for ( i=0; i < var->getNumberOfDirichletBC(); i++ ) {
            SolutionLib::FemDirichletBC *bc1 = var->getDirichletBC(i);
            bc1->setup(var->getCurrentOrder());
            std::vector<size_t> &list_bc_nodes = bc1->getListOfBCNodes();
            std::vector<double> &list_bc_values = bc1->getListOfBCValues();
            
	        // now loop over this vector
		    for ( j=0; j < list_bc_nodes.size(); j++ )
		    {
				size_t node_id = list_bc_nodes[j]; 
				double node_value = list_bc_values[j]; 
                std::map<size_t, ReductionKinNodeInfo*>::iterator my_bc; 
				my_bc = _bc_info.find( node_id ); 
                // now check whether this node has already been in the BC info
                if ( my_bc != _bc_info.end() )
                {   // this node is already there
					my_bc->second->set_comp_conc(i_var, node_value); 	             
                }
                else  // create a new structure and fill it in 
				{
					my_bc = _bc_info.begin(); 
					ReductionKinNodeInfo* bc_node = new ReductionKinNodeInfo( node_id, 
						                                                      this->_problem->getReductionScheme()->get_n_Comp(), 
						                                                      this->_problem->getReductionScheme()->get_n_eta_mob(), 
																			  this->_problem->getReductionScheme()->get_n_eta_immob(), 
																			  this->_problem->getReductionScheme()->get_n_xi_mob(), 
																			  this->_problem->getReductionScheme()->get_n_xi_immob(), 
																			  this->_problem->getReductionScheme() ); 
					bc_node->set_comp_conc( i_var, node_value ); 
					_bc_info.insert( my_bc, std::pair<size_t, ReductionKinNodeInfo*>(node_id, bc_node) ); 
				}  // end of if else
	        }  // end of for j
        }  // end of for i
    }  // end of for i_var


    // loop over all the boundary nodes, and 
	// transform these concentrations to eta and xi values
    std::map<size_t, ReductionKinNodeInfo*>::iterator bc_node_it; 
    std::vector<size_t> vec_bc_node_idx; 
    std::vector<std::vector<double>> vec_node_eta_values;
    std::vector<std::vector<double>> vec_node_xi_values; 

    for ( i=0; i < _linear_problem.size(); i++ )
    {
        std::vector<double> vec_eta_mob; 
        vec_node_eta_values.push_back(vec_eta_mob); 
    }
    
    for ( i=0; i < _non_linear_problem->getNumberOfVariables(); i++ )
    {
        std::vector<double> vec_xi;
        vec_node_xi_values.push_back(vec_xi); 
    }


    for ( bc_node_it = _bc_info.begin(); bc_node_it != _bc_info.end(); bc_node_it++ )
    {
        bc_node_it->second->transform(); 
        size_t node_idx = bc_node_it->second->get_node_id(); 
        vec_bc_node_idx.push_back(node_idx);

        for ( i=0; i < _linear_problem.size(); i++ )
        {
            double eta_mob_value = bc_node_it->second->get_eta_mob_value(i);
            vec_node_eta_values[i].push_back(eta_mob_value); 
        }

        for ( i=0; i < _non_linear_problem->getNumberOfVariables(); i++ )
        {
			double xi_value = bc_node_it->second->get_xi_mob_value(i);
            vec_node_xi_values[i].push_back(xi_value); 
        }
    }

    // imposing BC for eta
    for ( i=0; i < _linear_problem.size(); i++ )
     	_linear_problem[i]->getVariable(0)->addDirichletBC( new SolutionLib::FemDirichletBC( vec_bc_node_idx,  vec_node_eta_values[i] ) ); 

    // imposing BC for xi
	for ( i=0; i < _non_linear_problem->getNumberOfVariables(); i++ )
		_non_linear_problem->getVariable(i)->addDirichletBC( new SolutionLib::FemDirichletBC( vec_bc_node_idx,  vec_node_xi_values[i] ) ); 
};

template <
	class T_USER_FUNCTION_DATA, 
    class T_USER_FEM_PROBLEM,
    class T_USER_LINEAR_PROBLEM,
    class T_USER_LINEAR_SOLUTION,
    class T_USER_NON_LINEAR_PROBLEM, 
	class T_USER_NON_LINEAR_SOLUTION 
    >
int SingleStepKinReduction<T_USER_FUNCTION_DATA, T_USER_FEM_PROBLEM, T_USER_LINEAR_PROBLEM, T_USER_LINEAR_SOLUTION, T_USER_NON_LINEAR_PROBLEM, T_USER_NON_LINEAR_SOLUTION>
    ::solveTimeStep(const NumLib::TimeStep &t_n1)
{
	size_t i; 

	// solving linear problems one after the other
	for ( i=0; i < _lin_solutions.size(); i++)
	{
		// print solving information
		INFO("--Solving linear equation for eta_%d: ", i ); 
		_lin_solutions[i]->solveTimeStep( t_n1 );

		// if solution is accepted, 
		// if ( _lin_solutions[i]->accepted() )
		_function_data->set_eta_mob_node_values( i, _lin_solutions[i]->getCurrentSolution(0) ); 
	}
	// calcuate the reaction rates on each node
	_function_data->update_node_kin_reaction_rates(); 

	// solving the non-linear problem
	INFO("--Solving non-linear equations for xi:"); 
	_nlin_solution->solveTimeStep( t_n1 ); 

	// getting result
    // xi_mob
	for ( i=0; i < _nlin_solution->getProblem()->getNumberOfVariables(); i++)
	    _function_data->set_xi_mob_node_values( i, _nlin_solution->getCurrentSolution(i) ); 
    // xi_immob
    _function_data->update_xi_immob_node_values(); 

    //// st
    //std::vector<size_t> list_st_eqs_id;
    //std::vector<double> list_st_val;
    //for (size_t i_var=0; i_var<n_var; i_var++) {
    //    MyVariable* var = _problem->getVariable(i_var);
    //    for (size_t i=0; i<var->getNumberOfNeumannBC(); i++) {
    //        SolutionLib::IFemNeumannBC *bc2 = var->getNeumannBC(i);
    //        bc2->setup(var->getCurrentOrder());
    //        std::vector<size_t> &list_bc_nodes = bc2->getListOfBCNodes();
    //        std::vector<double> &list_bc_values = bc2->getListOfBCValues();
    //        DiscreteLib::convertToEqsValues(_dofManager, i_var, msh_id, list_bc_nodes, list_bc_values, list_st_eqs_id, list_st_val);
    //    }
    //}
    //(*_x_st) = .0;
    //for (size_t i=0; i<list_st_eqs_id.size(); i++) {
    //    (*_x_st)[list_st_eqs_id[i]] = list_st_val[i];
    //}

    // setup functions
    //_f_linear->reset(&t_n1, _x_n0);
    //_f_r->reset(&t_n1, _x_n0, _x_st);
    //_f_dx->reset(&t_n1, _x_n0);

    //// initial guess
    //*_x_n1_0 = *_x_n0;
    //for (size_t i=0; i<list_bc1_eqs_id.size(); i++) {
    //    (*_x_n1_0)[list_bc1_eqs_id[i]] = list_bc1_val[i];
    //}

    // solve
    //_f_nonlinear->solve(*_x_n1_0, *_x_n1);

    //// distribute solution vector to local vector for each variable
    //for (size_t i=0; i<n_var; i++) {
    //    //SolutionVector* vec_var = _problem->getVariable(i)->getIC()->getNodalValues();
    //    DiscreteLib::setLocalVector(_dofManager, i, msh_id, *_x_n1, *_vec_u_n1[i]->getDiscreteData());
    //}

    return 0;
}


template <
	class T_USER_FUNCTION_DATA, 
    class T_USER_FEM_PROBLEM,
    class T_USER_LINEAR_PROBLEM,
    class T_USER_LINEAR_SOLUTION,
    class T_USER_NON_LINEAR_PROBLEM, 
	class T_USER_NON_LINEAR_SOLUTION 
    >
bool SingleStepKinReduction<T_USER_FUNCTION_DATA, T_USER_FEM_PROBLEM, T_USER_LINEAR_PROBLEM, T_USER_LINEAR_SOLUTION, T_USER_NON_LINEAR_PROBLEM, T_USER_NON_LINEAR_SOLUTION>
    :: isBCNode(size_t node_idx)
{    
	if ( _bc_info.find( node_idx ) != _bc_info.end() )
		return true;  // is BC node
	else
		return false;  // not BC node
}

}

#endif  // end of ifndef
