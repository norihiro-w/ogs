/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file ITransientPartitionedAlgorithm.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <vector>
#include "NumLib/Coupling/ICoupledProblem.h"

namespace NumLib
{
class UnnamedParameterSet;
struct ParameterProblemMappingTable;

/**
 * \brief
 */
class ITransientPartitionedAlgorithm
{
public:
    virtual ~ITransientPartitionedAlgorithm() {};

    /// solve coupled problems
    /// @param subproblems    a list of subproblems
    /// @param vars           a container for shared variables
    /// @param mapping        mapping data between subproblems and shared variables
    virtual int solve(std::vector<ICoupledSystem*> &subproblems, UnnamedParameterSet &vars_t_n, UnnamedParameterSet &vars_t_n1, ParameterProblemMappingTable &mapping) = 0;
};


}
