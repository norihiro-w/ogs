/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "VecMatBuilder.h"

#include "MathLib/LinAlg/LinAlgBuilder.h"
#include "MathLib/LinAlg/MatrixOption.h"

namespace AssemblerLib
{

MathLib::IVector* VecMatBuilder::generateVector(MathLib::LinAlgLibType libType, const MeshComponentMap &dof_info)
{
    auto ghosts = dof_info.getGlobalIndicesOfGhosts();
    return MathLib::LinAlgBuilder::generateVector(libType, dof_info.nonGhostSize(), !dof_info.isParallel(), &ghosts);
}

MathLib::IMatrix* VecMatBuilder::generateMatrix(MathLib::LinAlgLibType libType, const MeshComponentMap &dof_info)
{
    MathLib::MatrixOption opt;
    opt.is_global_size = !dof_info.isParallel();
    return MathLib::LinAlgBuilder::generateMatrix(libType, dof_info.nonGhostSize(), &opt);
}

} // AssemblerLib

