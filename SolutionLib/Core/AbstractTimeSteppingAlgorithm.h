/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file AbstractTimeSteppingAlgorithm.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <vector>

#include "NumLib/TimeStepping/ITransientSystem.h"
#include "NumLib/TimeStepping/Algorithms/ITimeStepAlgorithm.h"


namespace SolutionLib
{

/**
 * \brief Abstract class for time-stepping method
 */
class AbstractTimeSteppingAlgorithm : public NumLib::ITransientSystem
{
public:
    /// @param tim Time step function
    AbstractTimeSteppingAlgorithm(NumLib::ITimeStepAlgorithm &tim) : _tim(&tim) {}

    /// destructor
    virtual ~AbstractTimeSteppingAlgorithm() {}

    /// get the time step function
    /// @return Time step function
    NumLib::ITimeStepAlgorithm* getTimeStepFunction() const {return _tim;}

    /// suggest the next time step
    /// @param time_currrent current time step object
    /// @return the next time
    double suggestNext(const NumLib::TimeStep &/*time_current*/)
    {
        _tim->next();
        return _tim->getTimeStep().current();
    }

    /// return if the next time step is the same as the given time
    /// @param time the time step object
    /// @bool true if this process is active with the given time
    bool isAwake(const NumLib::TimeStep &time)
    {
        return _tim->getTimeStep()==time;
    }

    ///
    virtual void accept(const NumLib::TimeStep &)
    {
        _tim->accepted();
    }


private:
    NumLib::ITimeStepAlgorithm* _tim;
};


}
