/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file OutputTimingBuilder.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <string>
#include <vector>

#include "IOutputTiming.h"

class OutputTimingBuilder
{
public:
    IOutputTiming* create(const std::string &name, size_t n, std::vector<double> *vec_time = 0);
};
