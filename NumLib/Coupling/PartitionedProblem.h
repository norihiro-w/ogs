/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file PartitionedProblem.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <vector>
#include <iostream>

#include "NumLib/IOSystem/UnnamedParameterSet.h"
#include "NumLib/IOSystem/NamedIOSystem.h"

#include "ICoupledProblem.h"
#include "AbstractPartitionedProblem.h"
#include "ParameterProblemMappingTable.h"
#include "Algorithm/PartitionedAlgorithm.h"

namespace NumLib
{

/**
 * \brief Partitioned problem
 */
class PartitionedProblem : public AbstractPartitionedProblem<ICoupledSystem>
{
public:
    typedef size_t InternalID;
    typedef size_t ExternalKey;

    /// 
    PartitionedProblem() : _algorithm(NULL)
    {
    }

    ///
    virtual ~PartitionedProblem()
    {
    };

    ///
    virtual void setAlgorithm(IPartitionedAlgorithm &algo)
    {
        _algorithm = &algo;
    }

    /// solve this system
    virtual int solve();

private:
    IPartitionedAlgorithm *_algorithm;
};



}
