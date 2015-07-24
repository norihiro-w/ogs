/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef VECMATBUILDER_H_
#define VECMATBUILDER_H_

#include "MathLib/LinAlg/LinAlgLibType.h"
#include "MeshComponentMap.h"

namespace MathLib
{
class IVector;
class IMatrix;
}

namespace AssemblerLib
{

class VecMatBuilder
{
public:
    static MathLib::IVector* generateVector(MathLib::LinAlgLibType libType, const MeshComponentMap &dof_info);
    static MathLib::IMatrix* generateMatrix(MathLib::LinAlgLibType libType, const MeshComponentMap &dof_info);
};

} // MathLib

#endif //VECMATBUILDER_H_

