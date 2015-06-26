/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file TransientCoupledSystem.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include "NumLib/Coupling/ICoupledProblem.h"
#include "NumLib/TimeStepping/TimeStep.h"
#include "NumLib/TimeStepping/ITransientSystem.h"

namespace NumLib
{

/**
 * \brief Interface class of transient coupling problems
 */
class ITransientCoupledSystem : public NumLib::ICoupledSystem, public ITransientSystem
{
public:
    int solve(const TimeStep &time) 
    {
        setCurrentTime(time);
        return solve();
    }
private:
    int solve()
    {
        return solveTimeStep(getCurrentTime());
    }
};

}
