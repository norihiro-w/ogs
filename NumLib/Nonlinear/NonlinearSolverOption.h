/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file NonlinearSolverOption.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <string>

#include "BaseLib/CodingTools.h"

namespace NumLib
{

struct NonlinerSolverOption
{
    enum SolverType
    {
        LINEAR,
        PICARD,
        NEWTON,
        INVALID
    };

    SolverType solver_type;
    double error_tolerance;
    long max_iteration;

    NonlinerSolverOption()
    {
        solver_type = LINEAR;
        error_tolerance = 1.e-6;
        max_iteration = 500;
    }

    SolverType getSolverType(const std::string &str)
    {
        RETURN_ENUM_IF_SAME_STRING(LINEAR, str);
        RETURN_ENUM_IF_SAME_STRING(PICARD, str);
        RETURN_ENUM_IF_SAME_STRING(NEWTON, str);

        return INVALID;
    }
};

} //end

