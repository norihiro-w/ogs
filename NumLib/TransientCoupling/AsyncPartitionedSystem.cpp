/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file AsyncPartitionedSystem.cpp
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#include "AsyncPartitionedSystem.h"

#include <limits>
#include <logog/include/logog.hpp>

namespace NumLib
{

double AsyncPartitionedSystem::suggestNext(const TimeStep &time_current)
{
    double t;
    t = std::numeric_limits<double>::max();
    for (size_t i=0; i<_list_subproblems.size(); i++) {
        ITransientCoupledSystem *solution = _list_subproblems[i];
        double suggest = solution->suggestNext(time_current);
        if (suggest>.0)
            t = std::min(t, suggest);
    }
    return t;
}

void AsyncPartitionedSystem::getActiveProblems(const TimeStep &time, std::vector<NumLib::ICoupledSystem*> &list_active_problems)
{
    for (size_t i=0; i<_list_subproblems.size(); i++) {
        ITransientCoupledSystem *solution = _list_subproblems[i];
        if (solution->isAwake(time)) {
            list_active_problems.push_back(solution);
        }
    }
}

int AsyncPartitionedSystem::solveTimeStep(const TimeStep &time)
{

    // copy previous time step result to current one
    //_vars_t_n.assign(_vars_t_n1);
    NumLib::UnnamedParameterSet *vars_t_n1 = getParameters();
    vars_t_n1->move(_vars_t_n);

    // list active problems
    std::vector<NumLib::ICoupledSystem*> list_active_problems;
    getActiveProblems(time, list_active_problems);

    if (list_active_problems.size()>0) {
        INFO("Solving a partitioned system with %d active sub problems...", list_active_problems.size());

        // solve
        for (size_t i=0; i<list_active_problems.size(); i++) {
            ITransientCoupledSystem *solution = static_cast<ITransientCoupledSystem*>(list_active_problems[i]);
            solution->setCurrentTime(time);
        }
        // set parameter state
        const size_t n_vars = vars_t_n1->size();
        std::vector<bool> list_org_state(n_vars);
        for (size_t i=0; i<n_vars; i++) {
            list_org_state[i] = vars_t_n1->isFixed(i);
            const NumLib::ParameterProblemMappingTable::PairSysVarId &v = _map._map_paraId2subproblem[i];
            const NumLib::ICoupledSystem *tmp_problem = v.first;
            bool is_active = (std::find(list_active_problems.begin(), list_active_problems.end(), tmp_problem)!=list_active_problems.end());
            if (!is_active) {
                vars_t_n1->setFixed(i, true);
            }
        }

        //if (list_active_problems.size()==1) {
        //    list_active_problems[0]->solve();
        //} else {
            _algorithm->solve(list_active_problems, _vars_t_n, *vars_t_n1, _map);
        //}

            // restore parameter state
            for (size_t i=0; i<n_vars; i++) {
                vars_t_n1->setFixed(i, list_org_state[i]);
            }

    }


    return 0;
};

bool AsyncPartitionedSystem::isAwake(const TimeStep &time)
{
    for (size_t i=0; i<_list_subproblems.size(); i++) {
        ITransientCoupledSystem *solution = _list_subproblems[i];
        if (solution->isAwake(time))
            return true;
    }
    return false;
};

void AsyncPartitionedSystem::accept(const TimeStep &time)
{
    for (size_t i=0; i<_list_subproblems.size(); i++) {
        ITransientCoupledSystem *solution = _list_subproblems[i];
        solution->accept(time);
    }
}

//int AsyncPartitionedSystem::find(const ITransientCoupledSystem &sub) const
//{
//    std::vector<ITransientCoupledSystem*>::const_iterator itr = std::find(_list_subproblems.begin(), _list_subproblems.end(), &sub);
//    if (itr!=_list_subproblems.end()) {
//        return (itr - _list_subproblems.begin());
//    } else {
//        return -1;
//    }
//}
//
///// check consistency
//bool AsyncPartitionedSystem::check() const
//{
//    bool flag = true;
//
//    // check for subproblems
//    for (size_t i=0; i<_list_subproblems.size(); i++) {
//        if (!_list_subproblems[i]->check())
//            flag = false;
//    }
//
//    // check consistency in input and output of subproblems
//    for (size_t i=0; i<_list_subproblems.size(); i++) {
//        const ITransientCoupledSystem* subproblem = _list_subproblems[i];
//        // check input parameters required for the subproblem
//        const NumLib::ParameterProblemMappingTable::ListOfInputVar &vec_registered_input_var = _map._list_subproblem_input_source[i];
//
//        std::vector<size_t> list1(vec_registered_input_var.size());
//        for (size_t j=0; j<list1.size(); j++) {
//            list1[j] = vec_registered_input_var[j].first;
//        }
//        std::vector<size_t> list2(subproblem->getNumberOfInputParameters());
//        for (size_t j=0; j<list2.size(); j++) {
//            list2[j] = j; //subproblem->getInputParameterKey(j);
//        }
//        sort(list1.begin(), list1.end());
//        sort(list2.begin(), list2.end());
//        std::vector<size_t> list_diff(list1.size()+list2.size());
//        std::vector<size_t>::iterator it = set_symmetric_difference(list1.begin(), list1.end(), list2.begin(), list2.end(), list_diff.begin());
//
//        int n_diff = it-list_diff.begin();
//        if (n_diff > 0) {
//            for (int k=0; k<n_diff; k++) {
//                std::cout << "*** Error: Inconsistent parameter found in subproblem " << i << " and param " << list_diff[k] << std::endl;
//            }
//
//            flag = false;
//        }
//    }
//
//    return flag;
//}
//
//size_t AsyncPartitionedSystem::addInputParameter(const std::string &name)
//{
//    size_t n = registerInputParameter(name);
//    _map._map_paraId2subproblem.push_back(std::make_pair<ITransientCoupledSystem*,size_t>(0, 0));
//    return n;
//}
//
//size_t AsyncPartitionedSystem::addOutputParameter(const std::string &name, ITransientCoupledSystem& sub_problem, size_t para_id_in_sub_problem)
//{
//    size_t var_id = registerOutputParameter(name);
//    // set reference
//    const NumLib::Parameter* v = sub_problem.getOutput(para_id_in_sub_problem);
//    if (v!=0)
//        setOutput(var_id, v);
//    // make a link between this and sub-problem variable
//    _map._map_paraId2subproblem.push_back(std::make_pair(&sub_problem, para_id_in_sub_problem));
//    // update a list of sub-problems
//    if (find(sub_problem)<0) {
//        addSubProblem(sub_problem);
//    }
//
//    return var_id;
//}
//
//void AsyncPartitionedSystem::connectInput(const std::string &this_para_name, ITransientCoupledSystem &subproblem, size_t subproblem_para_id)
//{
//    int this_para_id = SystemWithInOutParameters::getOutputParameterID(this_para_name);
//    if (this_para_id<0) return;
//
//    int subproblem_id = find(subproblem);
//    if (subproblem_id>=0) {
//        _map._list_subproblem_input_source[subproblem_id].push_back(std::make_pair(subproblem_para_id, this_para_id));
//    }
//}
//
//size_t AsyncPartitionedSystem::addSubProblem(ITransientCoupledSystem &sub_problem)
//{
//    _list_subproblems.push_back(&sub_problem);
//    _map._list_subproblem_input_source.resize(_list_subproblems.size());
//
//    return _list_subproblems.size();
//}
} //end
