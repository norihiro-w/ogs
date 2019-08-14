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

#include "NumLib/Fem/ShapeMatrixPolicy.h"

#include "ProcessLib/LIE/Common/HMatrixUtils.h"
#include "ProcessLib/LIE/ThermoHydroMechanics/ThermoHydroMechanicsProcessData.h"

#include "ThermoHydroMechanicsLocalAssemblerInterface.h"
#include "IntegrationPointDataFracture.h"

namespace ProcessLib
{
namespace LIE
{
namespace ThermoHydroMechanics
{
template <typename ShapeFunctionDisplacement,
          typename ShapeFunctionPressure,
          typename IntegrationMethod,
          int GlobalDim>
class ThermoHydroMechanicsLocalAssemblerFracture
    : public ThermoHydroMechanicsLocalAssemblerInterface
{
public:
    ThermoHydroMechanicsLocalAssemblerFracture(
        ThermoHydroMechanicsLocalAssemblerFracture const&) = delete;
    ThermoHydroMechanicsLocalAssemblerFracture(
        ThermoHydroMechanicsLocalAssemblerFracture&&) = delete;

    ThermoHydroMechanicsLocalAssemblerFracture(
        MeshLib::Element const& e,
        std::size_t const n_variables,
        std::size_t const local_matrix_size,
        std::vector<unsigned> const& dofIndex_to_localIndex,
        bool const is_axially_symmetric,
        unsigned const integration_order,
        ThermoHydroMechanicsProcessData<GlobalDim>& process_data);

    void preTimestepConcrete(std::vector<double> const& /*local_x*/,
                             double const /*t*/,
                             double const /*delta_t*/) override
    {
        for (auto& data : _ip_data)
        {
            data.pushBackState();
        }
    }

    void computeSecondaryVariableConcreteWithVector(
        const double t, Eigen::VectorXd const& local_x) override;

    Eigen::Map<const Eigen::RowVectorXd> getShapeMatrix(
        const unsigned integration_point) const override
    {
        auto const& N = _ip_data[integration_point].N_p;

        // assumes N is stored contiguously in memory
        return Eigen::Map<const Eigen::RowVectorXd>(N.data(), N.size());
    }

private:
    void assembleWithJacobianConcrete(double const t,
                                      Eigen::VectorXd const& local_x,
                                      Eigen::VectorXd const& local_xdot,
                                      Eigen::VectorXd& local_b,
                                      Eigen::MatrixXd& local_J) override;

    void assembleBlockMatricesWithJacobian(
        double const t,
        Eigen::Ref<const Eigen::VectorXd> const& p,
        Eigen::Ref<const Eigen::VectorXd> const& p_dot,
        Eigen::Ref<const Eigen::VectorXd> const& T,
        Eigen::Ref<const Eigen::VectorXd> const& T_dot,
        Eigen::Ref<const Eigen::VectorXd> const& g,
        Eigen::Ref<const Eigen::VectorXd> const& g_dot,
        Eigen::Ref<Eigen::VectorXd>
            rhs_p,
        Eigen::Ref<Eigen::VectorXd>
            rhs_T,
        Eigen::Ref<Eigen::VectorXd>
            rhs_g,
        Eigen::Ref<Eigen::MatrixXd>
            J_pp,
        Eigen::Ref<Eigen::MatrixXd>
            J_pT,
        Eigen::Ref<Eigen::MatrixXd>
            J_pg,
        Eigen::Ref<Eigen::MatrixXd>
            J_TT,
        Eigen::Ref<Eigen::MatrixXd>
            J_Tp,
        Eigen::Ref<Eigen::MatrixXd>
            J_Tg,
        Eigen::Ref<Eigen::MatrixXd>
            J_gg,
        Eigen::Ref<Eigen::MatrixXd>
            J_gp,
        Eigen::Ref<Eigen::MatrixXd>
            J_gT);

    // Types for displacement.
    using ShapeMatricesTypeDisplacement =
        ShapeMatrixPolicyType<ShapeFunctionDisplacement, GlobalDim>;
    using HMatricesType =
        HMatrixPolicyType<ShapeFunctionDisplacement, GlobalDim>;
    using HMatrixType = typename HMatricesType::HMatrixType;

    // Types for pressure.
    using ShapeMatricesTypePressure =
        ShapeMatrixPolicyType<ShapeFunctionPressure, GlobalDim>;

    // Types for the integration point data
    using IntegrationPointDataType =
        IntegrationPointDataFracture<HMatricesType,
                                     ShapeMatricesTypeDisplacement,
                                     ShapeMatricesTypePressure,
                                     GlobalDim>;

    ThermoHydroMechanicsProcessData<GlobalDim>& _process_data;
    std::vector<FractureProperty*> _fracture_props;
    std::vector<JunctionProperty*> _junction_props;
    std::unordered_map<int, int> _fracID_to_local;
    FractureProperty const* _fracture_property = nullptr;

    std::vector<IntegrationPointDataType,
                Eigen::aligned_allocator<IntegrationPointDataType>>
        _ip_data;

    Eigen::Vector3d _e_center_coords;

    static const int pressure_index_ = 0;
    static const int pressure_size_ = ShapeFunctionPressure::NPOINTS;
    static const int temperature_index_ = pressure_index_ + pressure_size_;
    static const int temperature_size_ = ShapeFunctionPressure::NPOINTS;
    static const int displacement_jump_index_ = temperature_index_ + temperature_size_;
    static const int displacement_jump_size_ =
        ShapeFunctionDisplacement::NPOINTS * GlobalDim;
};

}  // namespace ThermoHydroMechanics
}  // namespace LIE
}  // namespace ProcessLib

#include "ThermoHydroMechanicsLocalAssemblerFracture-impl.h"
