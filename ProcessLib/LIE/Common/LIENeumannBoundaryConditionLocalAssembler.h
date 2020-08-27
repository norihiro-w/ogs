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
class LIENeumannBoundaryConditionLocalAssembler final
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
    LIENeumannBoundaryConditionLocalAssembler(
        MeshLib::Element const& e,
        std::size_t const n_variables,
        std::size_t const local_matrix_size,
        std::vector<unsigned> const& dofIndex_to_localIndex,
        bool const is_axially_symmetric,
        unsigned const integration_order,
        ParameterLib::Parameter<double> const& neumann_bc_parameter,
        std::vector<std::unique_ptr<FractureProperty>> const& fracture_properties,
        std::vector<JunctionProperty> const& /*junction_properties*/,
        std::vector<unsigned> const& frac_ids)
        : Base(e, is_axially_symmetric, integration_order),
          _dofIndex_to_localIndex(dofIndex_to_localIndex),
          _local_rhs(n_variables * ShapeFunction::NPOINTS),
          //_local_rhs(n_variables * ShapeFunction::NPOINTS * GlobalDim),
          _local_rhs_active(local_matrix_size),
          _neumann_bc_parameter(neumann_bc_parameter),
          _element(e),
          _frac_ids(frac_ids)
    {
        for (auto fid : frac_ids)
        {
            _fracID_to_local.insert({fid, _fracture_props.size()});
            _fracture_props.push_back(&*fracture_properties[fid]);
        }

        // for (auto jid : vec_ele_connected_junctionIDs[e.getID()])
        // {
        //     _junction_props.push_back(&junction_properties[jid]);
        // }
    }

    void assemble(std::size_t const /*id*/,
                  NumLib::LocalToGlobalIndexMap const& dof_table_boundary,
                  double const t, const GlobalVector& /*x*/,
                  GlobalMatrix& /*K*/, GlobalVector& b,
                  GlobalMatrix* /*Jac*/) override
    {
        _local_rhs.setZero();

        unsigned const n_integration_points =
            Base::_integration_method.getNumberOfPoints();

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
            std::vector<double> const levelsets(uGlobalEnrichments(\
                _fracture_props, _junction_props, _fracID_to_local,
                pt));

            for (unsigned ifrac=0; ifrac<_fracID_to_local.size(); ifrac++)
            {
                _local_rhs.noalias() += ip_data.N * levelsets[ifrac] *
                                        parameter_node_values.dot(ip_data.N) *
                                        ip_data.weight;

            }
        }

        for (unsigned i = 0; i < _local_rhs_active.size(); i++)
        {
            _local_rhs_active[i] = _local_rhs[_dofIndex_to_localIndex[i]];
        }

        auto const indices = NumLib::getIndices(_element.getID(), dof_table_boundary);
        b.add(indices, _local_rhs_active);
    }

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

private:
    std::vector<unsigned> const _dofIndex_to_localIndex;
    Eigen::VectorXd _local_rhs;
    Eigen::VectorXd _local_rhs_active;
    ParameterLib::Parameter<double> const& _neumann_bc_parameter;
    MeshLib::Element const& _element;
    std::vector<unsigned> const _frac_ids;
    std::vector<FractureProperty*> _fracture_props;
    std::vector<JunctionProperty*> _junction_props;
    std::unordered_map<int, int> _fracID_to_local;

public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

}  // namespace LIE
}  // namespace ProcessLib
