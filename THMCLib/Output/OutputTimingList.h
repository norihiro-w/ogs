/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 */

#pragma once

#include <vector>
#include <algorithm>

#include "NumLib/TimeStepping/TimeStep.h"
#include "IOutputTiming.h"

class OutputTimingList : public IOutputTiming
{
public:
    explicit OutputTimingList(std::vector<double> &vec_time)
    : _list_time(vec_time)
    {
    }

    virtual ~OutputTimingList() {}

    virtual bool isActive(const NumLib::TimeStep &current_timestep)
    {
        return (_list_time.end() != std::find(_list_time.begin(), _list_time.end(), current_timestep.current()));
    }

private:
    std::vector<double> _list_time;
};
