/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file SerialStaggeredMethod.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include "AbstractPartitionedStaggeredMethod.h"
#include "BlockGaussSeidelMethod.h"

namespace NumLib
{

class SerialStaggeredMethod : public AbstractPartitionedStaggeredMethod<NumLib::BlockGaussSeidelMethod >
{
public:
    SerialStaggeredMethod() {};
    SerialStaggeredMethod(double epsilon, size_t max_count) : AbstractPartitionedStaggeredMethod<BlockGaussSeidelMethod >(epsilon, max_count) {};
    //SerialStaggeredMethod(IConvergenceCheck &checker, double epsilon, size_t max_count) : AbstractPartitionedStaggeredMethod<BlockGaussSeidelMethod >(checker, epsilon, max_count) {};
    virtual ~SerialStaggeredMethod() {};
};

}
