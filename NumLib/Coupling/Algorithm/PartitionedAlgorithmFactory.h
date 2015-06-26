/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file PartitionedAlgorithmFactory.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <string>
#include "BlockGaussSeidelMethod.h"
#include "BlockJacobiMethod.h"
#include "IConvergenceCheck.h"

namespace NumLib
{

class PartitionedAlgorithmFactory
{
public:
//    template<class T>
    static IPartitionedAlgorithm* create(const std::string &name, size_t max_itr, double epsilon)
    {
        if (name.compare("Jacobi")==0) {
            return new BlockJacobiMethod(epsilon, max_itr);
        } else if (name.compare("Gauss")==0) {
            return new BlockGaussSeidelMethod(epsilon, max_itr);
        }
        return 0;
    };
//    static IPartitionedAlgorithm* create(const std::string &name, IConvergenceCheck *checker, size_t max_itr, double epsilon)
//    {
//        if (name.compare("Jacobi")==0) {
//            return new BlockJacobiMethod(*checker, epsilon, max_itr);
//        } else if (name.compare("Gauss")==0) {
//            return new BlockGaussSeidelMethod(*checker, epsilon, max_itr);
//        }
//        return 0;
//    };
};


} //end
