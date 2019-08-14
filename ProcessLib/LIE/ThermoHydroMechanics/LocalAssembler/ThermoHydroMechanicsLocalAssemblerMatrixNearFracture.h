/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <unordered_map>
#include <vector>

#include "ProcessLib/LIE/Common/FractureProperty.h"
#include "ProcessLib/LIE/Common/JunctionProperty.h"

#include "ThermoHydroMechanicsLocalAssemblerMatrix.h"

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
class ThermoHydroMechanicsLocalAssemblerMatrixNearFracture
    : public ThermoHydroMechanicsLocalAssemblerMatrix<ShapeFunctionDisplacement,
                                                ShapeFunctionPressure,
                                                IntegrationMethod,
                                                GlobalDim>
{
    using Base = ThermoHydroMechanicsLocalAssemblerMatrix<ShapeFunctionDisplacement,
                                                    ShapeFunctionPressure,
                                                    IntegrationMethod,
                                                    GlobalDim>;

public:
    ThermoHydroMechanicsLocalAssemblerMatrixNearFracture(
        ThermoHydroMechanicsLocalAssemblerMatrixNearFracture const&) = delete;
    ThermoHydroMechanicsLocalAssemblerMatrixNearFracture(
        ThermoHydroMechanicsLocalAssemblerMatrixNearFracture&&) = delete;

    ThermoHydroMechanicsLocalAssemblerMatrixNearFracture(
        MeshLib::Element const& e,
        std::size_t const n_variables,
        std::size_t const local_matrix_size,
        std::vector<unsigned> const& dofIndex_to_localIndex,
        bool const is_axially_symmetric,
        unsigned const integration_order,
        ThermoHydroMechanicsProcessData<GlobalDim>& process_data);

private:
    void assembleWithJacobianConcrete(double const t,
                                      Eigen::VectorXd const& local_u,
                                      Eigen::VectorXd const& local_udot,
                                      Eigen::VectorXd& local_b,
                                      Eigen::MatrixXd& local_J) override;

    void preTimestepConcrete(std::vector<double> const& /*local_x*/,
                             double const /*t*/,
                             double const /*delta_t*/) override
    {
        for (auto& ip_data : _ip_data)
        {
            ip_data.pushBackState();
        }
    }

    void computeSecondaryVariableConcreteWithVector(
        double const t, Eigen::VectorXd const& local_x) override;

    Eigen::VectorXd compute_total_u(std::vector<double> const& levelsets,
                                    Eigen::VectorXd const& local_x) const;

    using Base::_element;
    using Base::_ip_data;
    using Base::_process_data;
    using Base::displacement_index;
    using Base::displacement_size;
    using Base::kelvin_vector_size;
    using Base::pressure_index;
    using Base::pressure_size;
    using Base::temperature_index;
    using Base::temperature_size;
    using typename Base::BMatricesType;
    using typename Base::ShapeMatricesTypeDisplacement;

    static const int displacement_jump_index =
        displacement_index + displacement_size;
    static const int displacement_jump_size = displacement_size;

    std::vector<FractureProperty*> _fracture_props;
    std::vector<JunctionProperty*> _junction_props;
    std::unordered_map<int, int> _fracID_to_local;
    unsigned _n_enrich_var;
    Eigen::Vector3d const _e_center_coords;
    std::vector<double> _ele_levelsets;
    bool _hasActiveLevelset = false;
};

}  // namespace ThermoHydroMechanics
}  // namespace LIE
}  // namespace ProcessLib

#include "ThermoHydroMechanicsLocalAssemblerMatrixNearFracture-impl.h"
