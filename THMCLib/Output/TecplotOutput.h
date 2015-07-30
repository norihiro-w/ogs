/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 */

#pragma once

#include <string>

#include "BaseLib/OrderedMap.h"
#include "NumLib/TimeStepping/TimeStep.h"

#include "IOutput.h"

class TecplotOutput : public IOutput
{
public:
    virtual void write(const NumLib::TimeStep &current_time,
            BaseLib::OrderedMap<std::string, OutputVariableInfo> &data);

private:
    bool _init = true;
};
