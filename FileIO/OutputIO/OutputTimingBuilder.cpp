/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file OutputTimingBuilder.cpp
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#include "OutputTimingBuilder.h"

#include "OutputTimingStepPeriodic.h"

IOutputTiming* OutputTimingBuilder::create(const std::string &name, size_t n, std::vector<double>* /*vec_time*/)
{
    if (name.compare("STEPS")==0) {
        return new OutputTimingStepPeriodic(n);
    }
    return NULL;
}
