/**
 * \file
 * \author Norihiro Watanabe
 * \date   2013-08-13
 * \brief
 *
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef GAUSSLEGENDRETRI_H_
#define GAUSSLEGENDRETRI_H_

#include <array>

namespace MathLib
{

/// Gauss-Legendre quadrature on triangles
///
/// \tparam ORDER   integration order.
template <unsigned ORDER>
struct GaussLegendreTri {
    static const unsigned Order = ORDER;
    static const unsigned NPoints = ORDER;
    static const std::array<std::array<double, 3>, NPoints> X;
    static const double W[NPoints];
};

template <>
struct GaussLegendreTri<2> {
    static const unsigned Order = 2;
    static const unsigned NPoints = 3;
    static const std::array<std::array<double, 3>, NPoints> X;
    static const double W[NPoints];
};

template <>
struct GaussLegendreTri<3> {
    static const unsigned Order = 3;
    static const unsigned NPoints = 4;
    static const std::array<std::array<double, 3>, NPoints> X;
    static const double W[NPoints];
};

}

#endif //GAUSSLEGENDRETRI_H_
