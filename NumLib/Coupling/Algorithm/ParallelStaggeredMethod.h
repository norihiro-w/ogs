/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file ParallelStaggeredMethod.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include "AbstractPartitionedStaggeredMethod.h"
#include "BlockJacobiMethod.h"

namespace NumLib
{

class ParallelStaggeredMethod : public AbstractPartitionedStaggeredMethod<NumLib::BlockJacobiMethod >
{
public:
    ParallelStaggeredMethod() {};
    ParallelStaggeredMethod(double epsilon, size_t max_count) : AbstractPartitionedStaggeredMethod<BlockJacobiMethod >(epsilon, max_count) {};
    //ParallelStaggeredMethod(IConvergenceCheck &checker, double epsilon, size_t max_count) : AbstractPartitionedStaggeredMethod<BlockJacobiMethod >(checker, epsilon, max_count) {};
    virtual ~ParallelStaggeredMethod() {};
};

}
