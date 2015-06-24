/**
 * \file
 * \author Norihiro Watanabe
 * \date   2012-02-17
 * \brief Helper macros.
 *
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef CODINGTOOLS_H
#define CODINGTOOLS_H

#include <cstdlib>

#define DISALLOW_COPY_AND_ASSIGN(TypeName) \
    TypeName(const TypeName&);   \
    TypeName &operator=(const TypeName&)

#define RETURN_ENUM_IF_SAME_STRING(TypeName,str) \
    if (str.compare(#TypeName)==0) return TypeName;

namespace BaseLib
{
const std::size_t index_npos = -1;
}

#endif
