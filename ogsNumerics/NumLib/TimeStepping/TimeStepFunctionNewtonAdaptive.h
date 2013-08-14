/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file TimeStepFunctionNewtonAdaptive.h
 *
 * Created on 2013-08-07 by Haibing Shao and Norihiro Watanabe
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
    /**
      * This function takes the log information and 
      * transfer it into the current TimeStepFunction class. 
      * In this case, only the number of Newton-Raphson iterations
      * is stored. 
      */
    virtual void updateLog(BaseLib::Options &log);

    /**
      * If the time if moving forward, then update the size of 
      * last time step. It will be used to evaluate the next time step. 
      */
    virtual void accept(double t)
    {
        if (t > getPrevious())
            _dt_pre = t - getPrevious();
        AbstractTimeStepFunction::accept(t);
    }

    /**
      * Overwritten suggestNext() function. 
      * In this function, the next time step size is evaluated based
      * on the size and number of Newton-Raphson iterations in the last
      * time step. 
      * For example, last time step size is 1.0 sec and it took 7 iterations, 
      * then we loop over the _iter_times_vector, find the corresponding multiplier
      * is 0.5. Then the next suggested time step size will be 1.0 * 0.5 = 0.5 sec. 
      */
    virtual double suggestNext(double /*t_current*/);

private:
    /**
      * this vector stores the number of iterations
      * to which the respective multiplier coefficient will be applied. 
      */
    std::vector<int> _iter_times_vector;
    
    /**
      * this vector stores the multiplier coefficients
      */
    std::vector<double> _multiplier_vector; 

    /**
      * the minimum allowed time step size
      */
    double _min_ts; 
    
    /**
      * the maximum allowed time step size
      */
    double _max_ts; 

    /**
      * the number of Newton-Raphson iterations
      * from the previous time step
      */
    int _iter_times;

    /**
      * previous time step size
      */
    double _dt_pre; 
};


}
