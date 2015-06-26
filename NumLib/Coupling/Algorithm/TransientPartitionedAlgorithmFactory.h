/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file TransientPartitionedAlgorithmFactory.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <string>
#include "TransientPartitionedAlgorithm.h"
#include "IConvergenceCheck.h"

namespace NumLib
{

class TransientPartitionedAlgorithmFactory
{
public:
    static ITransientPartitionedAlgorithm* create(const std::string &name, size_t max_itr, double epsilon)
    {
        if (name.compare("Parallel")==0) {
            return new ParallelStaggeredMethod(epsilon, max_itr);
        } else if (name.compare("Serial")==0) {
            return new SerialStaggeredMethod(epsilon, max_itr);
        }
        return 0;
    };
//    static ITransientPartitionedAlgorithm* create(const std::string &name, IConvergenceCheck *checker, size_t max_itr, double epsilon)
//    {
//        if (name.compare("Parallel")==0) {
//            return new ParallelStaggeredMethod(*checker, epsilon, max_itr);
//        } else if (name.compare("Serial")==0) {
//            return new SerialStaggeredMethod(*checker, epsilon, max_itr);
//        }
//        return 0;
//    };
};


} //end
