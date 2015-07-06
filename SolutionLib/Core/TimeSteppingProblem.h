/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file TimeSteppingProblem.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include "BaseLib/CodingTools.h"

namespace NumLib
{
class ITimeStepAlgorithm;
}

namespace SolutionLib
{

class TimeSteppingProblem
{
public:
    TimeSteppingProblem() : _tim(nullptr) {}

    virtual ~TimeSteppingProblem()
    {
//        delete _tim;
    }

    /// set  a time stepping function
    void setTimeSteppingFunction(NumLib::ITimeStepAlgorithm* f)
    {
        _tim = f;
    }

    /// get this time stepping function
    NumLib::ITimeStepAlgorithm* getTimeSteppingFunction() const {return _tim;}

private:
    DISALLOW_COPY_AND_ASSIGN(TimeSteppingProblem);

    NumLib::ITimeStepAlgorithm* _tim;
};

} //end
