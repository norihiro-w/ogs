/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file IConvergenceCheck.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include "NumLib/IOSystem/UnnamedParameterSet.h"

namespace NumLib
{

/**
 *
 */
class IConvergenceCheck
{
public:
    virtual ~IConvergenceCheck() {};

    virtual bool isConverged(UnnamedParameterSet& vars_prev, UnnamedParameterSet& vars_current, double eps, double &v_diff) = 0;
};

} //end
