/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file TimeSteppingControllerWithOutput.cpp
 *
 * Created on 2012-07-31 by Norihiro Watanabe
 */

#include "TimeSteppingControllerWithOutput.h"

namespace ogs6
{

void TimeSteppingControllerWithOutput::doSomethingAfterTimeStepAccepted(const NumLib::TimeStep &t_current) const
{
    _out_controller->outputData(t_current);
}

}
