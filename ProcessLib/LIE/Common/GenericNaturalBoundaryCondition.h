/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "MeshLib/MeshSubset.h"
#include "ProcessLib/BoundaryCondition/BoundaryCondition.h"

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
class GenericNaturalBoundaryConditionLocalAssemblerInterface;

namespace LIE
{
template <typename BoundaryConditionData,
          template <typename, typename, unsigned>
          class LocalAssemblerImplementation>
class GenericNaturalBoundaryCondition final : public BoundaryCondition
{
public:
    /// Create a boundary condition process from given config,
    /// DOF-table, and a mesh subset for a given variable and its component.
    /// A local DOF-table, a subset of the given one, is constructed.
    template <typename Data>
    GenericNaturalBoundaryCondition(
        unsigned const integration_order, unsigned const shapefunction_order,
        NumLib::LocalToGlobalIndexMap const& dof_table_bulk,
        int const variable_id, int const component_id,
        unsigned const global_dim, MeshLib::Mesh const& bc_mesh, Data&& data,
        std::vector<FractureProperty*> const& fracture_props,
        std::vector<JunctionProperty*> const& junction_props,
        std::unordered_map<int, int> const& fracID_to_local);

    /// Calls local assemblers which calculate their contributions to the global
    /// matrix and the right-hand-side.
    void applyNaturalBC(const double t, GlobalVector const& x, GlobalMatrix& K,
                        GlobalVector& b, GlobalMatrix* Jac) override;

private:
    /// Data used in the assembly of the specific boundary condition.
    BoundaryConditionData _data;

    /// A lower-dimensional mesh on which the boundary condition is defined.
    MeshLib::Mesh const& _bc_mesh;
//    /// Vector of lower-dimensional elements on which the boundary condition is
//    /// defined.
//    std::vector<MeshLib::Element*> _elements;

//    std::unique_ptr<MeshLib::MeshSubset const> _mesh_subset_all_nodes;

    /// Local dof table, a subset of the global one restricted to the
    /// participating number of _elements of the boundary condition.
    std::unique_ptr<NumLib::LocalToGlobalIndexMap> _dof_table_boundary;

    /// Integration order for integration over the lower-dimensional elements
    //unsigned const _integration_order;

    /// Local assemblers for each element of #_elements.
    std::vector<
        std::unique_ptr<GenericNaturalBoundaryConditionLocalAssemblerInterface>>
        _local_assemblers;
};

}  // namespace LIE
}  // namespace ProcessLib

#include "GenericNaturalBoundaryCondition-impl.h"