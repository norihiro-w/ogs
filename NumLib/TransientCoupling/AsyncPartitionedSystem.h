/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file AsyncPartitionedSystem.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <vector>

#include "NumLib/Coupling/ICoupledProblem.h"
#include "NumLib/Coupling/MonolithicProblem.h"
#include "NumLib/Coupling/PartitionedProblem.h"
#include "NumLib/Coupling/Algorithm/TransientPartitionedAlgorithm.h"
#include "TransientCoupledSystem.h"

namespace NumLib
{

/**
 * \brief Asynchronized partitioned problem
 */
class AsyncPartitionedSystem : public NumLib::AbstractPartitionedProblem<ITransientCoupledSystem>
{
public:
    typedef size_t InternalID;
    typedef size_t ExternalKey;

    AsyncPartitionedSystem() : _algorithm(0)
    {
    }

    virtual ~AsyncPartitionedSystem() {};

    virtual void setAlgorithm(NumLib::ITransientPartitionedAlgorithm &algo)
    {
        _algorithm = &algo;
    }

    double suggestNext(const TimeStep &time_current);
    
    int solveTimeStep(const TimeStep &time);
    
    bool isAwake(const TimeStep &time);

    void accept(const TimeStep &time);

    void getActiveProblems(const TimeStep &time, std::vector<NumLib::ICoupledSystem*> &list_active_problems);

private:
    NumLib::ITransientPartitionedAlgorithm *_algorithm;
    std::vector<TimeStep> _list_synchronize_time;
    
    NumLib::UnnamedParameterSet _vars_t_n;
};


}
