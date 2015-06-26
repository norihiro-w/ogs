/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file IPartitionedAlgorithm.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <vector>
#include "NumLib/Coupling/ICoupledProblem.h"
//#include "IConvergenceCheck.h"

namespace NumLib
{
class UnnamedParameterSet;
struct ParameterProblemMappingTable;

/**
 * \brief Interface class for partitioned algorithm
 */
class IPartitionedAlgorithm
{
public:

    virtual ~IPartitionedAlgorithm() {};

//    virtual void setConvergenceCheck(IConvergenceCheck &checker) = 0;

    /// solve coupled problems
    /// @param subproblems    a list of subproblems
    /// @param vars           a container for shared variables
    /// @param mapping        mapping data between subproblems and shared variables
    virtual int solve(const std::vector<ICoupledSystem*> &subproblems, UnnamedParameterSet &vars, const ParameterProblemMappingTable &mapping) = 0;
};

}
