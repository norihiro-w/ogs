/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file NewtonCotes.h
 *
 * Created on 2012-06-25 by Norihiro Watanabe
 */

#ifndef MATHLIB_NEWTONCOTES_H
#define MATHLIB_NEWTONCOTES_H

#include <cstddef>

namespace MathLib
{

/**
 * \brief Newton Cotes quadrature method
 *
 */
class NewtonCotes
{
    /// midpoint integration
    double middpoint(double (*fun)(double), double a, double b, std::size_t n);

    /// trapezoid integration
    double trapezoid(double (*fun)(double), double a, double b, std::size_t n);

    /// simpson integration
    double simpson(double (*fun)(double), double a, double b, std::size_t n);
};

}

#endif  // MATHLIB_NEWTONCOTES_H
