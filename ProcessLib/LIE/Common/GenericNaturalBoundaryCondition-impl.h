/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "MeshLib/MeshSearch/NodeSearch.h"  // for getUniqueNodes

#include "ProcessLib/BoundaryCondition/GenericNaturalBoundaryConditionLocalAssembler.h"
#include "ProcessLib/Utils/CreateLocalAssemblers.h"

namespace ProcessLib
{
namespace LIE
{
template <typename BoundaryConditionData,
          template <typename, typename, unsigned>
          class LocalAssemblerImplementation>
template <typename Data>
GenericNaturalBoundaryCondition<BoundaryConditionData,
                                LocalAssemblerImplementation>::
    GenericNaturalBoundaryCondition(
        unsigned const integration_order, unsigned const shapefunction_order,
        NumLib::LocalToGlobalIndexMap const& dof_table_bulk,
        int const variable_id, int const component_id,
        unsigned const global_dim, MeshLib::Mesh const& bc_mesh, Data&& data,
        std::vector<std::unique_ptr<FractureProperty>> const& fracture_props,
        std::vector<JunctionProperty> const& junction_props,
        std::vector<unsigned> const& frac_ids)
    : _data(std::forward<Data>(data)), _bc_mesh(bc_mesh)
{
    static_assert(std::is_same<typename std::decay<BoundaryConditionData>::type,
                               typename std::decay<Data>::type>::value,
                  "Type mismatch between declared and passed BC data.");

    // check basic data consistency
    if (variable_id >=
            static_cast<int>(dof_table_bulk.getNumberOfVariables()) ||
        component_id >=
            dof_table_bulk.getNumberOfVariableComponents(variable_id))
    {
        OGS_FATAL(
            "Variable id or component id too high. Actual values: (%d, %d), "
            "maximum values: (%d, %d).",
            variable_id, component_id, dof_table_bulk.getNumberOfVariables(),
            dof_table_bulk.getNumberOfVariableComponents(variable_id));
    }

    if (_bc_mesh.getDimension() + 1 != global_dim)
    {
        OGS_FATAL(
            "The dimension of the given boundary mesh (%d) is not by one lower "
            "than the bulk dimension (%d).",
            _bc_mesh.getDimension(), global_dim);
    }

    if (!_bc_mesh.getProperties().template existsPropertyVector<std::size_t>(
            "bulk_node_ids"))
    {
        OGS_FATAL(
            "The required bulk node ids map does not exist in the boundary "
            "mesh '%s'.",
            _bc_mesh.getName().c_str());
    }

    std::vector<MeshLib::Node*> const& bc_nodes = _bc_mesh.getNodes();
    DBUG("Found %d nodes for Natural BCs for the variable %d and component %d",
         bc_nodes.size(), variable_id, component_id);

    MeshLib::MeshSubset bc_mesh_subset(_bc_mesh, bc_nodes);

    // Create local DOF table from the BC mesh subset for the given variable and
    // component id.
    _dof_table_boundary.reset(dof_table_bulk.deriveBoundaryConstrainedMap(
        variable_id, {component_id}, std::move(bc_mesh_subset)));

    createLocalAssemblers<LocalAssemblerImplementation>(
        global_dim, _bc_mesh.getElements(), *_dof_table_boundary,
        shapefunction_order, _local_assemblers, _bc_mesh.isAxiallySymmetric(),
        integration_order, _data, fracture_props, junction_props, frac_ids,
        variable_id);
}

// template <typename BoundaryConditionData,
//           template <typename, typename, unsigned>
//           class LocalAssemblerImplementation>
// template <typename Data>
// GenericNaturalBoundaryCondition<BoundaryConditionData,
//                                 LocalAssemblerImplementation>::
//     GenericNaturalBoundaryCondition(
//         unsigned const integration_order, unsigned const shapefunction_order,
//         NumLib::LocalToGlobalIndexMap const& dof_table_bulk,
//         int const variable_id, int const component_id,
//         unsigned const global_dim, MeshLib::Mesh const& bc_mesh, Data&& data,
// 		FractureProperty const& fracture_prop)
//     : _data(std::forward<Data>(data)), _bc_mesh(bc_mesh)
// //      _elements(std::move(elements)),
// //      _integration_order(integration_order)
// {
//     assert(component_id < dof_table_bulk.getNumberOfComponents());

//     std::vector<MeshLib::Node*> nodes = MeshLib::getUniqueNodes(_elements);
//     DBUG("Found %d nodes for Natural BCs for the variable %d and component %d",
//          nodes.size(), variable_id, component_id);

//     auto const& mesh_subsets =
//         dof_table_bulk.getMeshSubsets(variable_id, component_id);

//     // TODO extend the node intersection to all parts of mesh_subsets, i.e.
//     // to each of the MeshSubset in the mesh_subsets.
//     _mesh_subset_all_nodes.reset(
//         mesh_subsets.getMeshSubset(0).getIntersectionByNodes(nodes));
//     MeshLib::MeshSubsets all_mesh_subsets{_mesh_subset_all_nodes.get()};

//     // Create local DOF table from intersected mesh subsets for the given
//     // variable and component ids.
//     _dof_table_boundary.reset(dof_table_bulk.deriveBoundaryConstrainedMap(
//         variable_id, {component_id}, std::move(all_mesh_subsets), _elements));

//     createLocalAssemblers<LocalAssemblerImplementation>(
//         global_dim, _elements, *_dof_table_boundary, shapefunction_order,
//         _local_assemblers, is_axially_symmetric, _integration_order, _data,
//         fracture_prop, variable_id);
// }

//template <typename BoundaryConditionData,
//          template <typename, typename, unsigned>
//          class LocalAssemblerImplementation>
//GenericNaturalBoundaryCondition<
//    BoundaryConditionData,
//    LocalAssemblerImplementation>::~GenericNaturalBoundaryCondition()
//{
//    for (auto e : _elements)
//        delete e;
//}

template <typename BoundaryConditionData,
          template <typename, typename, unsigned>
          class LocalAssemblerImplementation>
void GenericNaturalBoundaryCondition<
    BoundaryConditionData,
    LocalAssemblerImplementation>::applyNaturalBC(const double t,
                                                  const GlobalVector& x,
                                                  GlobalMatrix& K,
                                                  GlobalVector& b,
                                                  GlobalMatrix* Jac)
{
    GlobalExecutor::executeMemberOnDereferenced(
        &GenericNaturalBoundaryConditionLocalAssemblerInterface::assemble,
        _local_assemblers, *_dof_table_boundary, t, x, K, b, Jac);
}

}  // namespace LIE
}  // namespace ProcessLib
