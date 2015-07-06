/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 */

#pragma once

#include <string>
#include <logog/include/logog.hpp>

#include "Process.h"

namespace SolutionLib
{
class AbstractTimeSteppingAlgorithm;
}

namespace THMCLib
{

/**
 * \brief Implementation of Process (ITransientSystem) class for monolithic system
 *
 */
class AbstractTransientProcess : public Process
{
public:
    ///
    AbstractTransientProcess(const std::string &pcs_name, size_t n_in_parameters, size_t n_out_parameters)
    : Process(pcs_name, n_in_parameters, n_out_parameters)
    {}

    ///
    virtual ~AbstractTransientProcess() {}

    ///
    virtual int solveTimeStep(const NumLib::TimeStep &time);

    /// 
    virtual double suggestNext(const NumLib::TimeStep &time_current) ;

    ///
    virtual bool isAwake(const NumLib::TimeStep &time);

    ///
    virtual void accept(const NumLib::TimeStep &time);

protected:
    ///
    virtual SolutionLib::AbstractTimeSteppingAlgorithm* getSolution() = 0;
    ///
    virtual void initializeTimeStep(const NumLib::TimeStep &/*time*/) {};
    ///
    virtual void postSolutionAlgorithm(const NumLib::TimeStep &/*time*/) {};
    ///
    virtual void updateOutputParameter(const NumLib::TimeStep &time) = 0;
    ///
    virtual void postTimeStep(const NumLib::TimeStep &/*time*/) {};
    ///
    virtual void output(const NumLib::TimeStep &time) = 0;
};

} //end ProcessLib

