/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <vector>

#include "NumLib/DOF/DOFTableUtil.h"
#include "NumLib/Extrapolation/ExtrapolatableElement.h"

#include "LocalAssemblerInterface.h"


namespace ProcessLib
{
namespace SmallDeformationWithPT
{

struct SmallDeformationWithPTGlobalAssembler
{
public:
    void assembleResidual(
        const std::size_t mesh_item_id,
        SmallDeformationWithPTLocalAssemblerInterface& local_assembler,
        const NumLib::LocalToGlobalIndexMap& dof_table,
        const double t, GlobalVector &r)
    {
        auto const indices = NumLib::getIndices(mesh_item_id, dof_table);
        auto local_r = r.get(indices);

        //_local_r_data.clear();
        local_assembler.assembleResidual(t, local_r);
    }

private:
    std::vector<double> _local_r_data;
};

}  // namespace SmallDeformationWithPT
}  // namespace ProcessLib
