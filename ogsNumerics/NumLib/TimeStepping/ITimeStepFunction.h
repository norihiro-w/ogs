/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file ITimeStepFunction.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include "BaseLib/Options.h"

#include "TimeUnit.h"

namespace NumLib
{

class ITimeStepFunction
{
public:
    virtual double getBeginning() const = 0;
    virtual double getEnd() const = 0;
    virtual double getPrevious() const = 0;
    virtual double getNext(double t_current) = 0;
    virtual void accept(double t_current) = 0;
    virtual ITimeStepFunction* clone() = 0;
    virtual ~ITimeStepFunction() {};
    virtual void updateLog(BaseLib::Options &log) = 0;
};

}
