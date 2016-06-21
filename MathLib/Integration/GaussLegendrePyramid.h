/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef GAUSSLEGENDREPYRAMID_H_
#define GAUSSLEGENDREPYRAMID_H_

#include <array>
#include <mathlib_export.h>

namespace MathLib
{

/// Gauss-Legendre quadrature on pyramid
///
/// \tparam ORDER   integration order.
template <unsigned ORDER>
struct GaussLegendrePyramid {
	static MATHLIB_EXPORT const unsigned Order = ORDER;
	static MATHLIB_EXPORT const unsigned NPoints = ORDER;
	static MATHLIB_EXPORT const std::array<std::array<double, 3>, NPoints> X;
	static MATHLIB_EXPORT const double W[NPoints];
};

template <>
struct GaussLegendrePyramid<2> {
	static MATHLIB_EXPORT const unsigned Order = 2;
	static MATHLIB_EXPORT const unsigned NPoints = 5;
	static MATHLIB_EXPORT const std::array<std::array<double, 3>, NPoints> X;
	static MATHLIB_EXPORT const double W[NPoints];
};

template <>
struct GaussLegendrePyramid<3> {
	static MATHLIB_EXPORT const unsigned Order = 3;
	static MATHLIB_EXPORT const unsigned NPoints = 13;
	static MATHLIB_EXPORT const std::array<std::array<double, 3>, NPoints> X;
	static MATHLIB_EXPORT const double W[NPoints];
};

}

#endif //GAUSSLEGENDREPYRAMID_H_
