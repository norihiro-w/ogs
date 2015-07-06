/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file AbstractTimeIndependentProcess.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <string>
#include <limits>
#include "Process.h"

namespace THMCLib
{

/**
 *
 */
class AbstractTimeIndependentProcess
: public Process
{
public:
    AbstractTimeIndependentProcess(const std::string &pcs_name, size_t n_in_parameters, size_t n_out_parameters)
    : Process(pcs_name, n_in_parameters, n_out_parameters)
    {
    }

    virtual ~AbstractTimeIndependentProcess() {};

    double suggestNext(const NumLib::TimeStep &/*time_current*/) { return std::numeric_limits<double>::max(); }

    bool isAwake(const NumLib::TimeStep &/*time*/) { return true;  }

    virtual void accept(const NumLib::TimeStep &/*time*/) {};

};

}

