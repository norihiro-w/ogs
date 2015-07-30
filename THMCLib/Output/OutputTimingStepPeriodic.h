/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file OutputTimingStepPeriodic.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include "NumLib/TimeStepping/TimeStep.h"
#include "IOutputTiming.h"

class OutputTimingStepPeriodic : public IOutputTiming
{
public:
    explicit OutputTimingStepPeriodic(size_t step)
    : _steps(step)
    {
    }

    virtual ~OutputTimingStepPeriodic() {};

    virtual bool isActive(const NumLib::TimeStep &current_timestep)
    {
        return (current_timestep.steps() % _steps == 0);
    }

private:
    size_t _steps;
};
