/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "ProcessLib/BoundaryCondition/CreateBoundaryCondition.h"
#include "JunctionProperty.h"

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
}
}  // namespace ProcessLib

namespace ProcessLib
{
namespace LIE
{
/// A boundary condition builder for displacement jumps. Boundary
/// integration, e.g. for Neumann BC, should take into account the leveset
/// function.
class BoundaryConditionBuilder : public ProcessLib::BoundaryConditionBuilder
{
public:
    explicit BoundaryConditionBuilder(
        std::vector<FractureProperty*> const& fracture_props,
        std::vector<JunctionProperty*> const& junction_props,
        std::unordered_map<int, int> const& fracID_to_local)
        : _fracture_props(fracture_props),
        _junction_props(junction_props),
        _fracID_to_local(fracID_to_local)
    {
    }

    virtual std::unique_ptr<BoundaryCondition> createBoundaryCondition(
        const BoundaryConditionConfig& config,
        const NumLib::LocalToGlobalIndexMap& dof_table,
        const MeshLib::Mesh& mesh, const int variable_id,
        const unsigned integration_order, const unsigned shapefunction_order,
        const std::vector<std::unique_ptr<ParameterLib::ParameterBase>>& parameters,
        const Process& process) override;

    std::vector<FractureProperty*> const& _fracture_props;
    std::vector<JunctionProperty*> const& _junction_props;
    std::unordered_map<int, int> const& _fracID_to_local;
};

}  // namespace LIE
}  // namespace ProcessLib
