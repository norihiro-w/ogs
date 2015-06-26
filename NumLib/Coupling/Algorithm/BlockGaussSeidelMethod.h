/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file BlockGaussSeidelMethod.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include "AbstractIterativePartitionedMethod.h"

namespace NumLib
{

/**
 * \brief Block Gauss-Seidel iterative partitioned method
 */
class BlockGaussSeidelMethod : public AbstractIterativePartitionedMethod
{
public:
    BlockGaussSeidelMethod() {};
    ///
    BlockGaussSeidelMethod(double epsilon, size_t max_count) : AbstractIterativePartitionedMethod(epsilon, max_count)
    {
    }
    //BlockGaussSeidelMethod(IConvergenceCheck &checker, double epsilon, size_t max_count) : AbstractIterativePartitionedMethod(checker, epsilon, max_count)
//    {
//    }
    ///
    virtual ~BlockGaussSeidelMethod() {};

protected:
    bool isFixed() const {return true;};
    ///
    void doPostAfterSolve( ICoupledSystem & problem, UnnamedParameterSet& parameter_table, const ParameterProblemMappingTable &mapping ) //;
//    void BlockGaussSeidelMethod::doPostAfterSolve( ICoupledSystem & problem, UnnamedParameterSet& parameter_table, const ParameterProblemMappingTable &mapping )
    {
        // update shared variables
        AbstractIterativePartitionedMethod::updateParameterTable(problem, mapping, isFixed(), parameter_table);
    }
};



}
