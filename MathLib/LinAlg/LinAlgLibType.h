/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef LINALGLIBTYPE_H_
#define LINALGLIBTYPE_H_

namespace MathLib
{

enum class LinAlgLibType
{
    Dense,
    Eigen,
    Lis,
    EigenLis,
    PETSc
};

} // MathLib

#endif //LINALGLIBTYPE_H_

