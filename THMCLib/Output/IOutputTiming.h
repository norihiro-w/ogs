/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file IOutputTiming.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include "NumLib/TimeStepping/TimeStep.h"

class IOutputTiming
{
public:
    virtual ~IOutputTiming() {};
    virtual bool isActive(const NumLib::TimeStep &current_timestep) = 0;
};

