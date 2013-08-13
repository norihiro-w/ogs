/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file AbstractTimeStepFunction.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include "TimeUnit.h"
#include "ITimeStepFunction.h"

namespace NumLib
{

class AbstractTimeStepFunction : public ITimeStepFunction
{
public:
    AbstractTimeStepFunction(double t0, double tn, TimeUnit::type t_unit)
    {
        _t0 = TimeUnit::convertToSec(t_unit, t0);
        _tn = TimeUnit::convertToSec(t_unit, tn);
        _t_previous = _t_next = _t0;
        _steps = 0;
    };
    virtual ~AbstractTimeStepFunction() {};
    double getBeginning() const {return _t0;};
    double getEnd() const {return _tn;};
    void setEnd(double t) {_tn = t;};
    virtual void accept(double t_current)
    {
        if (t_current == _t_next && _t_previous < _t_next) {
            _t_previous = _t_next;
            _steps++;
        }
    }
    double getPrevious() const {return _t_previous;};
    double getNext(double t_current)
    {
        _t_next = suggestNext(t_current);
        return _t_next;
    }
    size_t getStep() const {return _steps;};

    virtual void updateLog(BaseLib::Options &/*log*/) {};

protected:
    virtual double suggestNext(double t_current) = 0;
private:
    double _t0;
    double _tn;
    double _t_previous;
    double _t_next;
    size_t _steps;
};


}
