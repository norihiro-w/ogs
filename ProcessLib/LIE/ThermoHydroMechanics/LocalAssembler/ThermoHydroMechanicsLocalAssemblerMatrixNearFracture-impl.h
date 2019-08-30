
#pragma once

#include "ThermoHydroMechanicsLocalAssemblerMatrixNearFracture.h"

namespace ProcessLib
{
namespace LIE
{
namespace ThermoHydroMechanics
{

namespace
{
bool hasActiveLevelset(std::vector<double> const& levelsets)
{
    for (auto const v : levelsets)
        if (v != 0)
            return true;
    return false;
}
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int GlobalDim>
ThermoHydroMechanicsLocalAssemblerMatrixNearFracture<ShapeFunctionDisplacement,
                                               ShapeFunctionPressure,
                                               IntegrationMethod, GlobalDim>::
    ThermoHydroMechanicsLocalAssemblerMatrixNearFracture(
        MeshLib::Element const& e,
        std::size_t const n_variables,
        std::size_t const local_matrix_size,
        std::vector<unsigned> const& dofIndex_to_localIndex,
        bool const is_axially_symmetric,
        unsigned const integration_order,
        ThermoHydroMechanicsProcessData<GlobalDim>& process_data)
    : ThermoHydroMechanicsLocalAssemblerMatrix<ShapeFunctionDisplacement,
                                         ShapeFunctionPressure,
                                         IntegrationMethod, GlobalDim>(
          e, n_variables, local_matrix_size, dofIndex_to_localIndex,
          is_axially_symmetric, integration_order, process_data),
         _e_center_coords(e.getCenterOfGravity().getCoords())
{
    for (auto fid : process_data.vec_ele_connected_fractureIDs[e.getID()])
    {
        _fracID_to_local.insert({fid, _fracture_props.size()});
        _fracture_props.push_back(&*_process_data.fracture_properties[fid]);
    }

    for (auto jid : process_data.vec_ele_connected_junctionIDs[e.getID()])
    {
        _junction_props.push_back(&_process_data.junction_properties[jid]);
    }
    _n_enrich_var = _fracture_props.size() + _junction_props.size();
    _ele_levelsets = uGlobalEnrichments(_fracture_props, _junction_props, _fracID_to_local, _e_center_coords);
    _hasActiveLevelset = hasActiveLevelset(_ele_levelsets);
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int GlobalDim>
void ThermoHydroMechanicsLocalAssemblerMatrixNearFracture<
    ShapeFunctionDisplacement, ShapeFunctionPressure, IntegrationMethod,
    GlobalDim>::assembleWithJacobianConcrete(double const t,
                                             Eigen::VectorXd const& local_x,
                                             Eigen::VectorXd const& local_x_dot,
                                             Eigen::VectorXd& local_b,
                                             Eigen::MatrixXd& local_J)
{
    auto p = const_cast<Eigen::VectorXd&>(local_x).segment(pressure_index,
                                                           pressure_size);
    auto p_dot = const_cast<Eigen::VectorXd&>(local_x_dot)
                     .segment(pressure_index, pressure_size);
    if (_process_data.deactivate_matrix_in_flow)
    {
        Base::setPressureOfInactiveNodes(t, p);
        Base::setPressureDotOfInactiveNodes(p_dot);
    }
    auto T = const_cast<Eigen::VectorXd&>(local_x).segment(temperature_index,
        temperature_size);
    auto T_dot = const_cast<Eigen::VectorXd&>(local_x_dot)
        .segment(temperature_index, temperature_size);
    auto const u = local_x.segment(displacement_index, displacement_size);
    auto const u_dot =
        local_x_dot.segment(displacement_index, displacement_size);

    auto rhs_p = local_b.segment(pressure_index, pressure_size);
    auto rhs_T = local_b.segment(temperature_index, temperature_size);
    auto rhs_u = local_b.segment(displacement_index, displacement_size);

    auto J_pp = local_J.block(pressure_index, pressure_index, pressure_size,
                              pressure_size);
    auto J_pT = local_J.block(pressure_index, temperature_index, pressure_size,
        temperature_size);
    auto J_pu = local_J.block(pressure_index, displacement_index, pressure_size,
                              displacement_size);
    auto J_Tp = local_J.block(temperature_index, pressure_index, temperature_size,
        pressure_size);
    auto J_TT = local_J.block(temperature_index, temperature_index, temperature_size,
        temperature_size);
    auto J_Tu = local_J.block(temperature_index, displacement_index, temperature_size,
        displacement_size);
    auto J_up = local_J.block(displacement_index, pressure_index,
                              displacement_size, pressure_size);
    auto J_uT = local_J.block(displacement_index, temperature_index,
        displacement_size, temperature_size);
    auto J_uu = local_J.block(displacement_index, displacement_index,
                              displacement_size, displacement_size);

    // levelset value of the element
    // remark: this assumes the levelset function is uniform within an element
    if (!_hasActiveLevelset)
    {
        // no DoF exists for displacement jumps. do the normal assembly
        Base::assembleBlockMatricesWithJacobian(t, p, p_dot, T, T_dot, u, u_dot, rhs_p, rhs_T,
                                                rhs_u, J_pp, J_pT, J_pu, J_TT, J_Tp, J_Tu, J_uu, J_up, J_uT);
        return;
    }

    // Displacement jumps should be taken into account
    // compute true displacements
    const Eigen::VectorXd total_u = compute_total_u(_ele_levelsets, local_x);
    const Eigen::VectorXd total_u_dot = compute_total_u(_ele_levelsets, local_x_dot);

    // evaluate residuals and Jacobians for pressure and displacements
    Base::assembleBlockMatricesWithJacobian(t, p, p_dot, T, T_dot, total_u, total_u_dot,
                                            rhs_p, rhs_T, rhs_u, J_pp, J_pT, J_pu, J_TT, J_Tp, J_Tu, J_uu, J_up, J_uT);

    // compute residuals and Jacobians for displacement jumps
    auto const n_enrich_var = _n_enrich_var;
    for (unsigned i = 0; i < n_enrich_var; i++)
    {
        auto g1_offset = displacement_jump_index + displacement_jump_size * i;
        auto rhs_g = local_b.segment(g1_offset, displacement_jump_size);
        auto J_pg = local_J.block(pressure_index, g1_offset,
                                pressure_size, displacement_jump_size);
        auto J_Tg = local_J.block(temperature_index, g1_offset,
                                temperature_size, displacement_jump_size);
        auto J_ug = local_J.block(displacement_index, g1_offset,
                                displacement_size, displacement_jump_size);
        auto J_gp = local_J.block(g1_offset, pressure_index,
                                displacement_jump_size, pressure_size);
        auto J_gT = local_J.block(g1_offset, temperature_index,
                                displacement_jump_size, temperature_size);
        auto J_gu = local_J.block(g1_offset, displacement_index,
                                displacement_jump_size, displacement_size);

        rhs_g = _ele_levelsets[i] * rhs_u;
        J_pg = _ele_levelsets[i] * J_pu;
        J_Tg = _ele_levelsets[i] * J_Tu;
        J_ug = _ele_levelsets[i] * J_uu;
        J_gp = _ele_levelsets[i] * J_up;
        J_gT = _ele_levelsets[i] * J_uT;
        J_gu = _ele_levelsets[i] * J_uu;

        for (unsigned j = 0; j < n_enrich_var; j++)
        {
            auto g2_offset = displacement_jump_index + displacement_jump_size * j;
            auto J_gg = local_J.block(g1_offset, g2_offset,
                                    displacement_jump_size, displacement_jump_size);
            J_gg = _ele_levelsets[i] * _ele_levelsets[j] * J_uu;
        }
    }

}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int GlobalDim>
void ThermoHydroMechanicsLocalAssemblerMatrixNearFracture<
    ShapeFunctionDisplacement, ShapeFunctionPressure, IntegrationMethod,
    GlobalDim>::
    computeSecondaryVariableConcreteWithVector(double const t,
                                               Eigen::VectorXd const& local_x)
{
    auto p = const_cast<Eigen::VectorXd&>(local_x).segment(pressure_index,
                                                           pressure_size);
    if (_process_data.deactivate_matrix_in_flow)
    {
        Base::setPressureOfInactiveNodes(t, p);
    }
    auto T = const_cast<Eigen::VectorXd&>(local_x).segment(temperature_index,
        temperature_size);
    auto u = local_x.segment(displacement_index, displacement_size);

    // levelset value of the element
    // remark: this assumes the levelset function is uniform within an element
    if (!_hasActiveLevelset)
    {
        // no DoF exists for displacement jumps. do the normal assembly
        Base::computeSecondaryVariableConcreteWithBlockVectors(t, p, T, u);
        return;
    }

    // Displacement jumps should be taken into account

    // compute true displacements
    const Eigen::VectorXd total_u = compute_total_u(_ele_levelsets, local_x);

    // evaluate residuals and Jacobians for pressure and displacements
    Base::computeSecondaryVariableConcreteWithBlockVectors(t, p, T, total_u);
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int GlobalDim>
Eigen::VectorXd ThermoHydroMechanicsLocalAssemblerMatrixNearFracture<
    ShapeFunctionDisplacement, ShapeFunctionPressure, IntegrationMethod,
    GlobalDim>::compute_total_u(std::vector<double> const& levelsets,
                                Eigen::VectorXd const& local_x) const
{
    auto const u = local_x.segment(displacement_index, displacement_size);
    Eigen::VectorXd total_u = u;
    auto const n_enrich_var = _n_enrich_var;
    for (unsigned i = 0; i < n_enrich_var; i++)
    {
        auto sub =
            const_cast<Eigen::VectorXd&>(local_x)
                .segment<displacement_jump_size>(displacement_jump_index +
                                                 displacement_jump_size * i);
        total_u += levelsets[i] * sub;
    }
    return total_u;
}

}  // namespace ThermoHydroMechanics
}  // namespace LIE
}  // namespace ProcessLib
