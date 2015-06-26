/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file TimeSteppingController.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <cstdlib>

namespace NumLib
{
class ITransientSystem;
class TimeStep;

/**
 * \brief Time stepping controller
 */    
class TimeSteppingController
{
public:
    ///
    TimeSteppingController(): _time_begin(.0), _root_subsystems(0) {}
    ///
    virtual ~TimeSteppingController() {}

    /// set the starting time
    void setBeginning(double time_begin) {_time_begin = time_begin;}
    /// return the starting time
    double getBeginning() const {return _time_begin;}

    /// set a transient system
    void setTransientSystem(ITransientSystem &sys) {_root_subsystems = &sys;}
    /// get a transient system
    ITransientSystem* getTransientSystem() const {return _root_subsystems;}

    /// solve systems until the given time
    std::size_t solve(double time_end);

protected:
    virtual void doSomethingAfterTimeStepAccepted(const TimeStep &/*t_current*/) const {};

private:
    double _time_begin;
    ITransientSystem* _root_subsystems;
};

}
