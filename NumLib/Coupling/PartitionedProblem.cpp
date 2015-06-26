/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file PartitionedProblem.cpp
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#include "PartitionedProblem.h"

namespace NumLib
{

int PartitionedProblem::solve()
{
    return _algorithm->solve(_list_subproblems, *getParameters(), _map);
}

} //end
