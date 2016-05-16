/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef OGS_MATHLIB_TYPES_H_
#define OGS_MATHLIB_TYPES_H_

namespace MathLib
{
#ifdef INDEX_TYPE_INT
typedef int GlobalIndexType;
#else
typedef std::size_t GlobalIndexType;
#endif

typedef double RealType;
}

#endif // OGS_MATHLIB_TYPES_H_
