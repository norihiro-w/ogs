/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "ProcessLib/BoundaryCondition/BoundaryConditionBuilder.h"

namespace MeshLib
{
class Mesh;
}

namespace NumLib
{
class LocalToGlobalIndexMap;
}

namespace ProcessLib
{
namespace LIE
{
struct FractureProperty;
struct JunctionProperty;
}
}  // namespace ProcessLib

namespace ProcessLib
{
namespace LIE
{
/// A boundary condition builder for displacement jumps. Boundary
/// integration, e.g. for Neumann BC, should take into account the leveset
/// function.
class LIEBoundaryConditionBuilder : public ProcessLib::BoundaryConditionBuilder
{
public:
    LIEBoundaryConditionBuilder(
        std::vector<std::unique_ptr<FractureProperty>> const& fracture_props,
        std::vector<JunctionProperty> const& junction_props,
        std::vector<unsigned> const& frac_ids)
        : _fracture_props(fracture_props),
        _junction_props(junction_props),
        _frac_ids(frac_ids)
    {
    }

    virtual std::unique_ptr<BoundaryCondition> createBoundaryCondition(
        const BoundaryConditionConfig& config,
        const NumLib::LocalToGlobalIndexMap& dof_table,
        const MeshLib::Mesh& mesh, const int variable_id,
        const unsigned integration_order, const unsigned shapefunction_order,
        const std::vector<std::unique_ptr<ParameterLib::ParameterBase>>& parameters,
        const Process& process) override;

    std::vector<std::unique_ptr<FractureProperty>> const& _fracture_props;
    std::vector<JunctionProperty> const& _junction_props;
    std::vector<unsigned> const _frac_ids;
};

}  // namespace LIE
}  // namespace ProcessLib
