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
#include "ProcessLib/Deformation/BMatrixPolicy.h"

#include "ProcessLib/LIE/ThermoHydroMechanics/ThermoHydroMechanicsProcessData.h"

#include "ThermoHydroMechanicsLocalAssemblerInterface.h"
#include "IntegrationPointDataMatrix.h"

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
class ThermoHydroMechanicsLocalAssemblerMatrix
    : public ThermoHydroMechanicsLocalAssemblerInterface
{
public:
    ThermoHydroMechanicsLocalAssemblerMatrix(
        ThermoHydroMechanicsLocalAssemblerMatrix const&) = delete;
    ThermoHydroMechanicsLocalAssemblerMatrix(ThermoHydroMechanicsLocalAssemblerMatrix&&) =
        delete;

    ThermoHydroMechanicsLocalAssemblerMatrix(
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

    Eigen::Map<const Eigen::RowVectorXd> getShapeMatrix(
        const unsigned integration_point) const override
    {
        auto const& N = _ip_data[integration_point].N_u;

        // assumes N is stored contiguously in memory
        return Eigen::Map<const Eigen::RowVectorXd>(N.data(), N.size());
    }

protected:
    void assembleWithJacobianConcrete(double const t,
                                      Eigen::VectorXd const& local_x,
                                      Eigen::VectorXd const& local_x_dot,
                                      Eigen::VectorXd& local_rhs,
                                      Eigen::MatrixXd& local_Jac) override;

    void assembleBlockMatricesWithJacobian(
        double const t,
        Eigen::Ref<const Eigen::VectorXd> const& p,
        Eigen::Ref<const Eigen::VectorXd> const& p_dot,
        Eigen::Ref<const Eigen::VectorXd> const& T,
        Eigen::Ref<const Eigen::VectorXd> const& T_dot,
        Eigen::Ref<const Eigen::VectorXd> const& u,
        Eigen::Ref<const Eigen::VectorXd> const& u_dot,
        Eigen::Ref<Eigen::VectorXd>
            rhs_p,
        Eigen::Ref<Eigen::VectorXd>
            rhs_T,
        Eigen::Ref<Eigen::VectorXd>
            rhs_u,
        Eigen::Ref<Eigen::MatrixXd>
            J_pp,
        Eigen::Ref<Eigen::MatrixXd>
            J_pT,
        Eigen::Ref<Eigen::MatrixXd>
            J_pu,
        Eigen::Ref<Eigen::MatrixXd>
            J_TT,
        Eigen::Ref<Eigen::MatrixXd>
            J_Tp,
        Eigen::Ref<Eigen::MatrixXd>
            J_Tu,
        Eigen::Ref<Eigen::MatrixXd>
            J_uu,
        Eigen::Ref<Eigen::MatrixXd>
            J_up,
        Eigen::Ref<Eigen::MatrixXd>
            J_uT);

    void computeSecondaryVariableConcreteWithVector(
        double const t, Eigen::VectorXd const& local_x) override;

    void computeSecondaryVariableConcreteWithBlockVectors(
        double const t,
        Eigen::Ref<const Eigen::VectorXd> const& p,
        Eigen::Ref<const Eigen::VectorXd> const& T,
        Eigen::Ref<const Eigen::VectorXd> const& u);

    void setPressureOfInactiveNodes(
        double const t, Eigen::Ref<Eigen::VectorXd> p);
    void setPressureDotOfInactiveNodes(Eigen::Ref<Eigen::VectorXd> p_dot);

    // Types for displacement.
    using ShapeMatricesTypeDisplacement =
        ShapeMatrixPolicyType<ShapeFunctionDisplacement, GlobalDim>;
    using BMatricesType =
        BMatrixPolicyType<ShapeFunctionDisplacement, GlobalDim>;

    // Types for pressure.
    using ShapeMatricesTypePressure =
        ShapeMatrixPolicyType<ShapeFunctionPressure, GlobalDim>;

    using IntegrationPointDataType =
        IntegrationPointDataMatrix<BMatricesType,
                                   ShapeMatricesTypeDisplacement,
                                   ShapeMatricesTypePressure,
                                   GlobalDim,
                                   ShapeFunctionDisplacement::NPOINTS>;

    ThermoHydroMechanicsProcessData<GlobalDim>& _process_data;

    std::vector<IntegrationPointDataType,
                Eigen::aligned_allocator<IntegrationPointDataType>>
        _ip_data;

    static constexpr int pressure_index = 0;
    static constexpr int pressure_size = ShapeFunctionPressure::NPOINTS;
    static constexpr int temperature_index = pressure_index + pressure_size;
    static constexpr int temperature_size = pressure_size;
    static constexpr int displacement_index = temperature_index + temperature_size;
    static constexpr int displacement_size =
        ShapeFunctionDisplacement::NPOINTS * GlobalDim;
    static const int kelvin_vector_size =
        MathLib::KelvinVector::KelvinVectorDimensions<GlobalDim>::value;
};

}  // namespace ThermoHydroMechanics
}  // namespace LIE
}  // namespace ProcessLib

#include "ThermoHydroMechanicsLocalAssemblerMatrix-impl.h"
