/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file AbstractIterativePartitionedMethod.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <cstddef>
#include <iostream>
#include <cmath>

#include "logog/include/logog.hpp"

#include "NumLib/IOSystem/UnnamedParameterSet.h"
#include "NumLib/Coupling/ICoupledProblem.h"
#include "NumLib/Coupling/ParameterProblemMappingTable.h"
#include "IPartitionedAlgorithm.h"
#include "IConvergenceCheck.h"

namespace NumLib
{

/**
 * \brief Abstract class for iterative partitioned methods
 */
//template <class T_CONVERGENCE_CHECK>
class AbstractIterativePartitionedMethod : public IPartitionedAlgorithm
{
public:
    ///
    AbstractIterativePartitionedMethod() : _max_itr(100), _itr_count(0), _epsilon(1e-3)/*, _convergence_check(0)*/ {}
    ///
    AbstractIterativePartitionedMethod(double epsilon, size_t max_count) : _max_itr(max_count), _itr_count(0), _epsilon(epsilon)/*, _convergence_check(0)*/ {}
    //AbstractIterativePartitionedMethod(IConvergenceCheck &checker, double epsilon, size_t max_count) : _max_itr(max_count), _itr_count(0), _epsilon(epsilon), _convergence_check(&checker) {}
    ///
    virtual ~AbstractIterativePartitionedMethod() {};

    ///
    //void setConvergenceCheck(IConvergenceCheck &checker) {_convergence_check = &checker;};

    /// get max. iteration
    size_t getMaximumIterationCounts() const {return _max_itr;};

    /// get epsilon
    double getEpsilon() const {return _epsilon;};

    /// get iteration count
    size_t getIterationCounts() const {return _itr_count;};

    /// solve
    //virtual int solve(const std::vector<ICoupledSystem*> &list_coupled_problems, UnnamedParameterSet &parameter_table, const ParameterProblemMappingTable &mapping);
    //template <class T_CONVERGENCE_CHECK>
    //int AbstractIterativePartitionedMethod<T_CONVERGENCE_CHECK>::solve (
    int solve (
            const std::vector<ICoupledSystem*> &list_coupled_problems,
            UnnamedParameterSet &parameter_table,
            const ParameterProblemMappingTable &mapping
            )
    {
//        if (_convergence_check==0) {
//            ERR("***Error: No convergence checker is specified in AbstractIterativePartitionedMethod::solve()");
//            return -1;
//        }
        const size_t n_subproblems = list_coupled_problems.size();

        // initialize variables
        updateParameterTable(mapping, isFixed(), parameter_table);

        // iteration start
        const size_t max_itr = getMaximumIterationCounts();
        bool is_converged = false;
        size_t i_itr = 0;
        double v_diff = .0;
        UnnamedParameterSet prev_parameter_table;
        do {
            //prev_parameter_table.assign(parameter_table);
            parameter_table.move(prev_parameter_table);

            // compute each solution
            size_t count_calculated = 0;
            for (size_t i=0; i<n_subproblems; i++) {
                ICoupledSystem *problem = list_coupled_problems[i];
                const std::vector<ParameterProblemMappingTable::PairInputVar> &problem_parameters = mapping._list_subproblem_input_source[i];

                // calculate all anyway in the 1st iteration
                //if (i_itr>0 && !isInputParametersUpdated(parameter_table, problem_parameters, problem))
                //    continue;

                setInputParameters(parameter_table, problem_parameters, problem);
                problem->solve();
                doPostAfterSolve(*problem, parameter_table, mapping);
                count_calculated++;
            }

            if (count_calculated>0) {
                doPostAfterSolveAll(parameter_table, mapping);
                parameter_table.finishUpdate();
            }

            if (n_subproblems>1) {
                // check convergence
                is_converged = true;
                for (size_t i=0; i<n_subproblems; i++) {
                    ICoupledSystem *problem = list_coupled_problems[i];
                    IConvergenceCheck* checker = problem->getConvergenceChecker();
                    double temp_diff = .0;
                    bool is_this_converged = checker->isConverged(prev_parameter_table, parameter_table, getEpsilon(), temp_diff);
                    v_diff = std::max(temp_diff, v_diff);
                    if (!is_this_converged) is_converged = false;
                }
//                if (_convergence_check!=0) {
//                    //T_CONVERGENCE_CHECK check;
//                    //is_converged = check.isConverged(prev_parameter_table, parameter_table, getEpsilon(), v_diff);
//                    is_converged = _convergence_check->isConverged(prev_parameter_table, parameter_table, getEpsilon(), v_diff);
//                }
                //std::cout.setf(std::ios_base::scientific, std::ios_base::floatfield);
                //std::cout.precision(3);
                //std::cout << v_diff << " ";
            } else {
                is_converged = true;
            }
            ++i_itr;
        } while (!is_converged && i_itr<max_itr);

        _itr_count = i_itr;

        if (n_subproblems>1) {
            INFO("------------------------------------------------------------------");
            INFO("*** Partitioned coupling solver computation result");
            INFO("nr. of subproblems : %d", n_subproblems);
            if (max_itr==1) {
                INFO("status             : iteration not required");
            } else {
                INFO("status             : %s", (is_converged ? "converged" : "***ERROR - DID NOT CONVERGE!"));
            }
            INFO("iteration          : %d/%d", i_itr, max_itr);
            INFO("residuals          : %1.3e (tolerance=%1.3e)", v_diff, getEpsilon());
            INFO("------------------------------------------------------------------");
        }

        //if (i_itr==max_itr) {
        //    std::cout << "the iteration reached the maximum count " << max_itr << std::endl;
        //}

        return 0;
    }
protected:
    virtual void doPostAfterSolve( ICoupledSystem& /*solution*/, UnnamedParameterSet& /*vars*/, const ParameterProblemMappingTable& /*mapping*/ )  {}

    virtual void doPostAfterSolveAll( UnnamedParameterSet& /*vars*/, const ParameterProblemMappingTable& /*mapping*/ ) {}

    virtual bool isFixed() const = 0;

    /// update parameter table from all problems
    void updateParameterTable(const ParameterProblemMappingTable &mapping, bool /*is_fixed*/, UnnamedParameterSet &parameter_table)
    {
        const size_t n_parameters = parameter_table.size();
        for (size_t i=0; i<n_parameters; i++) {
            const ParameterProblemMappingTable::PairSysVarId &sysvarid = mapping._map_paraId2subproblem[i];
            const ICoupledSystem *solution = sysvarid.first;
            if (solution!=0) {
                //if (is_fixed)
                //    parameter_table.setFixed(i, true);
                parameter_table.set(i, solution->getOutput(sysvarid.second));
            }
        }
    }

    /// update parameter table from one problem
    void updateParameterTable(const ICoupledSystem &src_problem, const ParameterProblemMappingTable &mapping, bool is_fixed, UnnamedParameterSet &parameter_table)
    {
        const size_t n_vars = parameter_table.size();
        for (size_t i=0; i<n_vars; i++) {
            const ParameterProblemMappingTable::PairSysVarId &v = mapping._map_paraId2subproblem[i];
            const ICoupledSystem *tmp_problem = v.first;
            if (tmp_problem==&src_problem) {
                if (is_fixed)
                    parameter_table.setFixed(i, true);
                parameter_table.set(i, src_problem.getOutput(v.second));
            }
        }
    }

    /// update problem parameter from parameter table
    void setInputParameters( const UnnamedParameterSet &parameter_table, const std::vector<ParameterProblemMappingTable::PairInputVar> &problem_parameters, ICoupledSystem* problem)
    {
        const size_t n_input_parameters = problem_parameters.size();
        for (size_t j=0; j<n_input_parameters; j++) {
            size_t local_var_id = problem_parameters[j].first;
            size_t shared_var_id = problem_parameters[j].second;
            problem->setInput(local_var_id, parameter_table.get(shared_var_id));
        }
    }

    /// is para
    bool isInputParametersUpdated( const UnnamedParameterSet &parameter_table, const std::vector<ParameterProblemMappingTable::PairInputVar> &problem_parameters, ICoupledSystem* /*problem*/)
    {
        const size_t n_input_parameters = problem_parameters.size();
        for (size_t j=0; j<n_input_parameters; j++) {
            //size_t local_var_id = problem_parameters[j].first;
            size_t shared_var_id = problem_parameters[j].second;
            if (parameter_table.isUpdated(shared_var_id)) return true;
        }
        return false;
    }

private:
    size_t _max_itr;
    size_t _itr_count;
    double _epsilon;
    //IConvergenceCheck *_convergence_check;
};




} //end
