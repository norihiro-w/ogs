/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <memory>

namespace BaseLib
{
class ConfigTree;
}

namespace NumLib
{
class TimeStepAlgorithm;

std::unique_ptr<TimeStepAlgorithm> createExternalTimeSteppingController(
    BaseLib::ConfigTree const& config);
}  // end of namespace NumLib
