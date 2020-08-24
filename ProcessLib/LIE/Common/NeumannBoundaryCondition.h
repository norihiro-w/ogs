/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "ProcessLib/BoundaryCondition/BoundaryCondition.h"
#include "GenericNaturalBoundaryCondition.h"
#include "NeumannBoundaryConditionLocalAssembler.h"
#include "ParameterLib/Parameter.h"

namespace MeshLib
{
class Element;
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
using NeumannBoundaryCondition =
    GenericNaturalBoundaryCondition<ParameterLib::Parameter<double> const&,
                                    NeumannBoundaryConditionLocalAssembler>;

std::unique_ptr<NeumannBoundaryCondition> createNeumannBoundaryCondition(
    BaseLib::ConfigTree const& config, MeshLib::Mesh const& bc_mesh,
    NumLib::LocalToGlobalIndexMap const& dof_table, int const variable_id,
    int const component_id, unsigned const integration_order,
    unsigned const shapefunction_order, unsigned const global_dim,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const&
        parameters,
    std::vector<FractureProperty*> const& fracture_props,
    std::vector<JunctionProperty*> const& junction_props,
    std::unordered_map<int, int> const& fracID_to_local);

// std::unique_ptr<BoundaryCondition> createNeumannBoundaryCondition(
//     BaseLib::ConfigTree const& config, MeshLib::Mesh const& bc_mesh,
//     std::vector<MeshLib::Element*>&& elements,
//     NumLib::LocalToGlobalIndexMap const& dof_table, int const variable_id,
//     int const component_id, bool is_axially_symmetric,
//     unsigned const integration_order, unsigned const shapefunction_order,
//     unsigned const global_dim,
//     std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
//     FractureProperty const& fracture_prop);

}  // namespace LIE
}  // namespace ProcessLib
