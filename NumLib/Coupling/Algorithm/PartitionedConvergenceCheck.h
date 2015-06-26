/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file PartitionedConvergenceCheck.h
 *
 * Created on 2012-08-17 by Norihiro Watanabe
 */

#pragma once

#include <vector>

#include "NumLib/Coupling/ParameterProblemMappingTable.h"
#include "IConvergenceCheck.h"

namespace NumLib
{

/**
 *
 */
template <class T_PROBLEM>
class PartitionedConvergenceCheck : public IConvergenceCheck
{
public:
    PartitionedConvergenceCheck(std::vector<T_PROBLEM*>* list_subproblems, ParameterProblemMappingTable* mapt)
    : _list_subproblems(list_subproblems), _map(mapt)
    {};

    virtual ~PartitionedConvergenceCheck() {};

    virtual bool isConverged(UnnamedParameterSet& vars_prev, UnnamedParameterSet& vars_current, double eps, double &v_diff)
    {
        bool is_converged = true;
        for (size_t i=0; i<_list_subproblems->size(); i++) {
            T_PROBLEM* prob = (*_list_subproblems)[i];
            IConvergenceCheck* check = prob->getConvergenceChecker();
            double temp_diff = .0;
            bool is_this_converged = check->isConverged(vars_prev, vars_current, eps, temp_diff);
            v_diff = std::max(temp_diff, v_diff);
            if (!is_this_converged) {
                is_converged = false;
            }
        }
        return is_converged;
    }

private:
    std::vector<T_PROBLEM*>* _list_subproblems;
    ParameterProblemMappingTable* _map;
};

} //end
