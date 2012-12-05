/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file GaussLegendre.h
 *
 * Created on 2012-06-25 by Norihiro Watanabe
 */

#ifndef MATHLIB_GAUSSLEGENDRE_H
#define MATHLIB_GAUSSLEGENDRE_H

#include <vector>
#include <cstddef>

namespace MathLib
{

/**
 * \brief Gauss-Ledgendre quadrature method
 *
 */
class GaussLegendre
{
public:
    /**
     * return a sampling point
     *
     * \param n_sample_points   the number of sampling points
     * \param point_id          point index
     * \return x
     */
    static double getPoint(std::size_t n_sample_points, std::size_t point_id);
    

    /**
     * return a weigth of the sampling point
     *
     * \param n_sample_points   the number of sampling points
     * \param point_id          point index
     * \return weight
     */
    static double getWeight(std::size_t n_sample_points, std::size_t point_id);


    /**
     * integrate the given function
     *
     * \param fun               integrant with one argument
     * \param sampling_size     the number of sampling points
     * \return integral value
     */
    double integrate(double (*fun)(double), int sampling_size);
};


}

#endif
