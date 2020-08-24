/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "NumLib/DOF/DOFTableUtil.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "ProcessLib/BoundaryCondition/GenericNaturalBoundaryConditionLocalAssembler.h"
#include "ParameterLib/Parameter.h"
#include "ProcessLib/Utils/InitShapeMatrices.h"
#include "ProcessLib/LIE/Common/LevelSetFunction.h"
#include "ProcessLib/LIE/Common/Utils.h"

namespace ProcessLib
{
namespace LIE
{
template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
class NeumannBoundaryConditionLocalAssembler final
    : public GenericNaturalBoundaryConditionLocalAssembler<
          ShapeFunction, IntegrationMethod, GlobalDim>
{
    using Base = GenericNaturalBoundaryConditionLocalAssembler<
        ShapeFunction, IntegrationMethod, GlobalDim>;
    using ShapeMatricesType = ShapeMatrixPolicyType<ShapeFunction, GlobalDim>;
    using NodalVectorType = typename Base::NodalVectorType;

public:
    /// The neumann_bc_value factor is directly integrated into the local
    /// element matrix.
    NeumannBoundaryConditionLocalAssembler(
        MeshLib::Element const& e,
        std::size_t const local_matrix_size,
        bool const is_axially_symmetric,
        unsigned const integration_order,
        ParameterLib::Parameter<double> const& neumann_bc_parameter,
        std::vector<FractureProperty*> const& fracture_props,
        std::vector<JunctionProperty*> const& junction_props,
        std::unordered_map<int, int> const& fracID_to_local,
        int variable_id)
        : Base(e, is_axially_symmetric, integration_order),
          _neumann_bc_parameter(neumann_bc_parameter),
          _local_rhs(local_matrix_size),
          _element(e),
          _fracture_props(fracture_props),
          _junction_props(junction_props),
          _fracID_to_local(fracID_to_local),
          _variable_id(variable_id)
    {
    }

    void assemble(std::size_t const id,
                  NumLib::LocalToGlobalIndexMap const& dof_table_boundary,
                  double const t, const GlobalVector& /*x*/,
                  GlobalMatrix& /*K*/, GlobalVector& b,
                  GlobalMatrix* /*Jac*/) override
    {
        _local_rhs.setZero();

        unsigned const n_integration_points =
            Base::_integration_method.getNumberOfPoints();

//        SpatialPosition pos;
//        pos.setElementID(id);

        // Get element nodes for the interpolation from nodes to integration
        // point.
        NodalVectorType parameter_node_values =
            _neumann_bc_parameter.getNodalValuesOnElement(Base::_element, t)
                .template topRows<ShapeFunction::MeshElement::n_all_nodes>();

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            auto const& ip_data = Base::_ns_and_weights[ip];

            // levelset functions
            auto const ip_physical_coords =
                 computePhysicalCoordinates(_element, ip_data.N);
            Eigen::Vector3d const pt(ip_physical_coords.getCoords());
            // double const levelsets = calculateLevelSetFunction(
            //      _fracture_prop, ip_physical_coords.getCoords());
            // std::vector<double> const levelsets(uGlobalEnrichments(
            //     e_fracture_props, e_junction_props, e_fracID_to_local, pt));
            std::vector<double> const levelsets(uGlobalEnrichments(\
                _fracture_props, _junction_props, _fracID_to_local, 
                pt));
            int fracid;

            _local_rhs.noalias() += ip_data.N * levelsets[fracid] *
                                    parameter_node_values.dot(ip_data.N) *
                                    ip_data.weight;
			
            // pos.setIntegrationPoint(ip);
            // auto const& sm = Base::_shape_matrices[ip];
            // auto const& wp = Base::_integration_method.getWeightedPoint(ip);

            // // levelset functions
            // auto const ip_physical_coords =
            //     computePhysicalCoordinates(_element, sm.N);
            // double const levelsets = calculateLevelSetFunction(
            //     _fracture_prop, ip_physical_coords.getCoords());

            // _local_rhs.noalias() += sm.N * levelsets *
            //                         _neumann_bc_parameter(t, pos)[0] * sm.detJ *
            //                         wp.getWeight() * sm.integralMeasure;
        }

        auto const indices = NumLib::getIndices(id, dof_table_boundary);
        b.add(indices, _local_rhs);
    }

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

private:
    ParameterLib::Parameter<double> const& _neumann_bc_parameter;
    NodalVectorType _local_rhs;
    MeshLib::Element const& _element;
    std::vector<FractureProperty*> const& _fracture_props;
    std::vector<JunctionProperty*> const& _junction_props;
    std::unordered_map<int, int> const& _fracID_to_local;
    int const _variable_id;

public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

}  // namespace LIE
}  // namespace ProcessLib
