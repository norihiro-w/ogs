/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 */

#include "OutputTimingBuilder.h"

#include "OutputTimingStepPeriodic.h"
#include "OutputTimingList.h"

IOutputTiming* OutputTimingBuilder::create(const std::string &name, size_t n, std::vector<double>* vec_time)
{
    if (name == "STEPS") {
        return new OutputTimingStepPeriodic(n);
    } else if (vec_time != nullptr) { //if (name == "LIST") {
        return new OutputTimingList(*vec_time);
    }
    return NULL;
}
