/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 */

#pragma once

namespace NumLib
{

/// Polynomial order
enum class PolynomialOrder
{
    Linear = 1,
    Quadratic = 2,
    Cubic = 3,
    Quartic = 4,
    INVALID = -1
};

} //end
