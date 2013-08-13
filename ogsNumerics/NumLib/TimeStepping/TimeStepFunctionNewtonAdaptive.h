/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file TimeStepFunctionNewtonAdaptive.h
 *
 * Created on 2013-08-07 by Haibing Shao
 */

#pragma once

#include <vector>

#include "TimeUnit.h"
#include "AbstractTimeStepFunction.h"

namespace NumLib
{


class TimeStepFunctionNewtonAdaptive : public AbstractTimeStepFunction
{
public:
	TimeStepFunctionNewtonAdaptive(double t0,
                                   double tn, 
                                   double min_ts, 
                                   double max_ts, 
                                   const std::vector<int> &iter_times_vector, 
                                   const std::vector<double> &multiplier_vector, 
                                   TimeUnit::type t_unit=TimeUnit::Second);

    TimeStepFunctionNewtonAdaptive* clone();

protected:
    virtual void updateLog(BaseLib::Options &log);

    virtual void accept(double t)
    {
        _dt_pre = t - getPrevious(); 
        AbstractTimeStepFunction::accept(t);
    }

    virtual double suggestNext(double /*t_current*/);


private:
    std::vector<int> _iter_times_vector;
    std::vector<double> _multiplier_vector; 
    double _min_ts; 
    double _max_ts; 
    int _iter_times;
    double _dt_pre;    // previous time step size
};


}
