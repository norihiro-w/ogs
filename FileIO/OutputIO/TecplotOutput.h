/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file TecplotOutput.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <string>
#include <vector>
#include <map>

#include "BaseLib/CodingTools.h"
#include "NumLib/Function/ITXFunction.h"
#include "NumLib/TimeStepping/TimeStep.h"

#include "IOutput.h"

class TecplotOutput : public IOutput
{
public:
    virtual void write(const NumLib::TimeStep &/*current_time*/, BaseLib::OrderedMap<std::string, OutputVariableInfo> &/*data*/)
    {

    }
};
