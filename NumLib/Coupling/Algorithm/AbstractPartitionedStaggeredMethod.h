/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file AbstractPartitionedStaggeredMethod.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <vector>
#include "NumLib/IOSystem/UnnamedParameterSet.h"
#include "NumLib/Coupling/ParameterProblemMappingTable.h"
#include "ITransientPartitionedAlgorithm.h"
#include "IConvergenceCheck.h"

namespace NumLib
{

/**
 * \brief
 */
template <class T_PARTITIONED>
class AbstractPartitionedStaggeredMethod : public ITransientPartitionedAlgorithm
{
public:
    AbstractPartitionedStaggeredMethod() : _part_method() {};
    AbstractPartitionedStaggeredMethod(double epsilon, size_t max_count) : _part_method(epsilon, max_count) {};
    //AbstractPartitionedStaggeredMethod(IConvergenceCheck &checker, double epsilon, size_t max_count) : _part_method(checker, epsilon, max_count) {};
    virtual ~AbstractPartitionedStaggeredMethod() {};

    /// solve
    int solve(std::vector<ICoupledSystem*> &subproblems, UnnamedParameterSet &parameters_t_n, UnnamedParameterSet &parameters_t_n1, ParameterProblemMappingTable &mapping)
    {
        initializeParameters(parameters_t_n, subproblems, mapping, parameters_t_n1);

        _part_method.solve(subproblems, parameters_t_n1, mapping);

        return 0;
    }

private:
    void initializeParameters(const UnnamedParameterSet &parameters_t_n, const std::vector<ICoupledSystem*> &/*subproblems*/, const ParameterProblemMappingTable &mapping, UnnamedParameterSet &parameters_t_n1)
    {
        const size_t n_para = parameters_t_n.size();
        for (size_t i=0; i<n_para; i++) {
            const ParameterProblemMappingTable::PairSysVarId &sysvarid = mapping._map_paraId2subproblem[i];
            const ICoupledSystem *solution = sysvarid.first;
            //size_t internalId = sysvarid.second;
            if (solution!=0)
                parameters_t_n1.set(i, parameters_t_n.get(i));
        }
    }
private:
    T_PARTITIONED _part_method;
};

}
