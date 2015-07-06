/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 */

#include "AbstractTransientProcess.h"

#include <logog/include/logog.hpp>

#include "NumLib/TransientCoupling/TransientMonolithicSystem.h"
#include "SolutionLib/Core/AbstractTimeSteppingAlgorithm.h"

namespace THMCLib
{

int AbstractTransientProcess::solveTimeStep(const NumLib::TimeStep &time)
{
    INFO("Solving %s...", getProcessName().c_str());
    initializeTimeStep(time);
    getSolution()->solveTimeStep(time);
    postSolutionAlgorithm(time);
    updateOutputParameter(time);
    return 0;
}

double AbstractTransientProcess::suggestNext(const NumLib::TimeStep &time_current)
{
    return getSolution()->suggestNext(time_current);
}

bool AbstractTransientProcess::isAwake(const NumLib::TimeStep &time)
{
    return getSolution()->isAwake(time);
}

void AbstractTransientProcess::accept(const NumLib::TimeStep &time)
{
    postTimeStep(time);
    output(time);
    getSolution()->accept(time);
}

}

