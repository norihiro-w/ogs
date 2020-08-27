/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <functional>
#include <memory>
#include <type_traits>
#include <typeindex>
#include <typeinfo>
#include <unordered_map>
#include <iostream>

#include "MeshLib/Elements/Elements.h"
#include "NumLib/DOF/LocalToGlobalIndexMap.h"
#include "NumLib/Fem/Integration/GaussLegendreIntegrationPolicy.h"

#ifndef OGS_MAX_ELEMENT_DIM
static_assert(false, "The macro OGS_MAX_ELEMENT_DIM is undefined.");
#endif

#ifndef OGS_MAX_ELEMENT_ORDER
static_assert(false, "The macro OGS_MAX_ELEMENT_ORDER is undefined.");
#endif

// The following macros decide which element types will be compiled, i.e.
// which element types will be available for use in simulations.

#ifdef OGS_ENABLE_ELEMENT_SIMPLEX
#define ENABLED_ELEMENT_TYPE_SIMPLEX 1u
#else
#define ENABLED_ELEMENT_TYPE_SIMPLEX 0u
#endif

#ifdef OGS_ENABLE_ELEMENT_CUBOID
#define ENABLED_ELEMENT_TYPE_CUBOID 1u << 1
#else
#define ENABLED_ELEMENT_TYPE_CUBOID 0u
#endif

#ifdef OGS_ENABLE_ELEMENT_PRISM
#define ENABLED_ELEMENT_TYPE_PRISM 1u << 2
#else
#define ENABLED_ELEMENT_TYPE_PRISM 0u
#endif

#ifdef OGS_ENABLE_ELEMENT_PYRAMID
#define ENABLED_ELEMENT_TYPE_PYRAMID 1u << 3
#else
#define ENABLED_ELEMENT_TYPE_PYRAMID 0u
#endif

// Dependent element types.
// Faces of tets, pyramids and prisms are triangles
#define ENABLED_ELEMENT_TYPE_TRI                                       \
    ((ENABLED_ELEMENT_TYPE_SIMPLEX) | (ENABLED_ELEMENT_TYPE_PYRAMID) | \
     (ENABLED_ELEMENT_TYPE_PRISM))
// Faces of hexes, pyramids and prisms are quads
#define ENABLED_ELEMENT_TYPE_QUAD                                     \
    ((ENABLED_ELEMENT_TYPE_CUBOID) | (ENABLED_ELEMENT_TYPE_PYRAMID) | \
     (ENABLED_ELEMENT_TYPE_PRISM))

// All enabled element types
#define OGS_ENABLED_ELEMENTS                                          \
    ((ENABLED_ELEMENT_TYPE_SIMPLEX) | (ENABLED_ELEMENT_TYPE_CUBOID) | \
     (ENABLED_ELEMENT_TYPE_PYRAMID) | (ENABLED_ELEMENT_TYPE_PRISM))

// Include only what is needed (Well, the conditions are not sharp).
#if OGS_ENABLED_ELEMENTS != 0
#include "NumLib/Fem/ShapeFunction/ShapeLine2.h"
#include "NumLib/Fem/ShapeFunction/ShapeLine3.h"
#endif

#if (OGS_ENABLED_ELEMENTS & ENABLED_ELEMENT_TYPE_TRI) != 0
#include "NumLib/Fem/ShapeFunction/ShapeTri3.h"
#include "NumLib/Fem/ShapeFunction/ShapeTri6.h"
#endif

#if (OGS_ENABLED_ELEMENTS & ENABLED_ELEMENT_TYPE_QUAD) != 0
#include "NumLib/Fem/ShapeFunction/ShapeQuad4.h"
#include "NumLib/Fem/ShapeFunction/ShapeQuad8.h"
#include "NumLib/Fem/ShapeFunction/ShapeQuad9.h"
#endif


namespace ProcessLib
{
namespace LIE
{
/// The LocalDataInitializer is a functor creating a local assembler data with
/// corresponding to the mesh element type shape functions and calling
/// initialization of the new local assembler data.
/// For example for MeshLib::Quad a local assembler data with template argument
/// NumLib::ShapeQuad4 is created.
template <typename LocalAssemblerInterface,
          template <typename, typename, int> class LocalAssemblerData,
          int GlobalDim, typename... ConstructorArgs>
class BCLocalDataInitializer final
{
public:
    using LADataIntfPtr = std::unique_ptr<LocalAssemblerInterface>;

    explicit BCLocalDataInitializer(
        NumLib::LocalToGlobalIndexMap const& dof_table,
        int const variable_id, int const component_id)
        : _dof_table(dof_table), _variable_id(variable_id), _component_id(component_id)
    {
        // REMARKS: At the moment, only a 2D mesh with 1D elements are
        // supported.

        // /// Quads and Hexahedra ///////////////////////////////////

#if (OGS_ENABLED_ELEMENTS & ENABLED_ELEMENT_TYPE_QUAD) != 0 && \
    OGS_MAX_ELEMENT_DIM >= 2 && OGS_MAX_ELEMENT_ORDER >= 1
        _builder[std::type_index(typeid(MeshLib::Quad))] =
            makeLocalAssemblerBuilder<NumLib::ShapeQuad4>();
#endif

#if (OGS_ENABLED_ELEMENTS & ENABLED_ELEMENT_TYPE_QUAD) != 0 && \
    OGS_MAX_ELEMENT_DIM >= 2 && OGS_MAX_ELEMENT_ORDER >= 2
        _builder[std::type_index(typeid(MeshLib::Quad8))] =
            makeLocalAssemblerBuilder<NumLib::ShapeQuad8>();
        _builder[std::type_index(typeid(MeshLib::Quad9))] =
            makeLocalAssemblerBuilder<NumLib::ShapeQuad9>();
#endif

        // /// Simplices ////////////////////////////////////////////////

#if (OGS_ENABLED_ELEMENTS & ENABLED_ELEMENT_TYPE_TRI) != 0 && \
    OGS_MAX_ELEMENT_DIM >= 2 && OGS_MAX_ELEMENT_ORDER >= 1
        _builder[std::type_index(typeid(MeshLib::Tri))] =
            makeLocalAssemblerBuilder<NumLib::ShapeTri3>();
#endif

#if (OGS_ENABLED_ELEMENTS & ENABLED_ELEMENT_TYPE_TRI) != 0 && \
    OGS_MAX_ELEMENT_DIM >= 2 && OGS_MAX_ELEMENT_ORDER >= 2
        _builder[std::type_index(typeid(MeshLib::Tri6))] =
            makeLocalAssemblerBuilder<NumLib::ShapeTri6>();
#endif

        // /// Lines ///////////////////////////////////

#if OGS_MAX_ELEMENT_DIM >= 2 && OGS_MAX_ELEMENT_ORDER >= 1
        _builder[std::type_index(typeid(MeshLib::Line))] =
            makeLocalAssemblerBuilder<NumLib::ShapeLine2>();
#endif

#if OGS_MAX_ELEMENT_DIM >= 2 && OGS_MAX_ELEMENT_ORDER >= 2
        _builder[std::type_index(typeid(MeshLib::Line3))] =
            makeLocalAssemblerBuilder<NumLib::ShapeLine3>();
#endif
    }

    /// Sets the provided \c data_ptr to the newly created local assembler data.
    ///
    /// \attention
    /// The index \c id is not necessarily the mesh item's id. Especially when
    /// having multiple meshes it will differ from the latter.
    void operator()(std::size_t const /*id*/,
                    MeshLib::Element const& mesh_item,
                    LADataIntfPtr& data_ptr,
                    ConstructorArgs&&... args) const
    {
        auto const type_idx = std::type_index(typeid(mesh_item));
        auto const it = _builder.find(type_idx);
        auto const eid = mesh_item.getID();

        if (it == _builder.end())
        {
            OGS_FATAL(
                "You are trying to build a local assembler for an unknown mesh "
                "element type (%s)."
                " Maybe you have disabled this mesh element type in your build "
                "configuration or this process requires higher order elements.",
                type_idx.name());
        }

        auto const n_local_dof = _dof_table.getNumberOfElementDOF(eid);
        //std::cout << _dof_table << std::endl;
        // auto const n_global_components =
        //     _dof_table.getNumberOfElementComponents(id);
        //auto const varIDs = _dof_table.getElementVariableIDs(eid);


        // INFO("n_local_dof = %d", n_local_dof);
        std::vector<unsigned> dofIndex_to_localIndex(n_local_dof, 0);
        //dofIndex_to_localIndex.resize(n_local_dof);
        unsigned dof_id = 0;
        unsigned local_id = 0;

        auto const& ms = _dof_table.getMeshSubset(0, 0);
//        auto const& ms = _dof_table.getMeshSubset(_variable_id, _component_id);
        auto const mesh_id = ms.getMeshID();

        for (unsigned k = 0; k < mesh_item.getNumberOfNodes(); k++)
        {
            MeshLib::Location l(mesh_id,
                                MeshLib::MeshItemType::Node,
                                mesh_item.getNodeIndex(k));
            auto global_index = _dof_table.getGlobalIndex(l, _variable_id, _component_id);
            if (global_index != NumLib::MeshComponentMap::nop)
            {
                dofIndex_to_localIndex[dof_id++] = local_id;
                // INFO("e=%d, k=%d, node=%d", mesh_item.getID(), k, mesh_item.getNodeIndex(k));
            }
            local_id++;
        }
        const unsigned nvar = 1; //TODO
        // INFO("n_local_dof = %d", n_local_dof);
        data_ptr = it->second(mesh_item, nvar, n_local_dof,
                              dofIndex_to_localIndex,
                              std::forward<ConstructorArgs>(args)...);
    }

private:
    using LADataBuilder = std::function<LADataIntfPtr(
        MeshLib::Element const& e,
        std::size_t const n_variables,
        std::size_t const local_matrix_size,
        std::vector<unsigned> const& dofIndex_to_localIndex,
        ConstructorArgs&&...)>;

    template <typename ShapeFunction>
    using IntegrationMethod = typename NumLib::GaussLegendreIntegrationPolicy<
        typename ShapeFunction::MeshElement>::IntegrationMethod;

    template <typename ShapeFunction>
    using LAData =
        LocalAssemblerData<ShapeFunction,
                                 IntegrationMethod<ShapeFunction>, GlobalDim>;

    /// A helper forwarding to the correct version of makeLocalAssemblerBuilder
    /// depending whether the global dimension is less than the shape function's
    /// dimension or not.
    template <typename ShapeFunction>
    static LADataBuilder makeLocalAssemblerBuilder()
    {
        return makeLocalAssemblerBuilder<ShapeFunction>(
            static_cast<std::integral_constant<
                bool, (GlobalDim >= ShapeFunction::DIM)>*>(nullptr));
    }

    /// Mapping of element types to local assembler constructors.
    std::unordered_map<std::type_index, LADataBuilder> _builder;

    NumLib::LocalToGlobalIndexMap const& _dof_table;
    int const _variable_id;
    int const _component_id;

    // local assembler builder implementations.
private:
    /// Generates a function that creates a new LocalAssembler of type
    /// LAData<ShapeFunction>. Only functions with shape function's dimension
    /// less or equal to the global dimension are instantiated, e.g.  following
    /// combinations of shape functions and global dimensions: (Line2, 1),
    /// (Line2, 2), (Line2, 3), (Hex20, 3) but not (Hex20, 2) or (Hex20, 1).
    template <typename ShapeFunction>
    static LADataBuilder makeLocalAssemblerBuilder(std::true_type*)
    {
        return [](MeshLib::Element const& e,
                  std::size_t const n_variables,
                  std::size_t const local_matrix_size,
                  std::vector<unsigned> const& dofIndex_to_localIndex,
                  ConstructorArgs&&... args)
        {
            //INFO("makeLocalAssemblerBuilder::local_matrix_size = %d", local_matrix_size);
            return LADataIntfPtr{new LAData<ShapeFunction>{
                e, n_variables, local_matrix_size, dofIndex_to_localIndex,
                std::forward<ConstructorArgs>(args)...}};
        };
    }

    /// Returns nullptr for shape functions whose dimensions are less than the
    /// global dimension.
    template <typename ShapeFunction>
    static LADataBuilder makeLocalAssemblerBuilder(std::false_type*)
    {
        return nullptr;
    }
};

}  // namespace LIE
}  // namespace ProcessLib

#undef ENABLED_ELEMENT_TYPE_SIMPLEX
#undef ENABLED_ELEMENT_TYPE_CUBOID
#undef ENABLED_ELEMENT_TYPE_PYRAMID
#undef ENABLED_ELEMENT_TYPE_PRISM
#undef ENABLED_ELEMENT_TYPE_TRI
#undef ENABLED_ELEMENT_TYPE_QUAD
#undef OGS_ENABLED_ELEMENTS
