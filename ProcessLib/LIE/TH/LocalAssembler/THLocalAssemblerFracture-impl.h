
#pragma once

#include <iostream>

#include "THLocalAssemblerFracture.h"

#include "MaterialLib/FractureModels/FractureIdentity2.h"

#include "ProcessLib/Utils/InitShapeMatrices.h"

#include "ProcessLib/LIE/Common/LevelSetFunction.h"

namespace ProcessLib
{
namespace LIE
{
namespace TH
{
template <int GlobalDim, typename RotationMatrix>
Eigen::Matrix<double, GlobalDim, GlobalDim> createRotatedTensor(
    RotationMatrix const& R, double const value)
{
    using M = Eigen::Matrix<double, GlobalDim, GlobalDim>;
    M tensor = M::Zero();
    tensor.diagonal().head(GlobalDim - 1).setConstant(value);
    return (R.transpose() * tensor * R).eval();
}

template <typename ShapeFunctionPressure, typename IntegrationMethod,
          int GlobalDim>
THLocalAssemblerFracture<ShapeFunctionPressure, IntegrationMethod, GlobalDim>::
    THLocalAssemblerFracture(
        MeshLib::Element const& e,
        std::size_t const /*n_variables*/,
        std::size_t const /*local_matrix_size*/,
        std::vector<unsigned> const& dofIndex_to_localIndex,
        bool const is_axially_symmetric,
        unsigned const integration_order,
        THProcessData<GlobalDim>& process_data)
    : THLocalAssemblerInterface(
          e, is_axially_symmetric,
          ShapeFunctionPressure::NPOINTS * 2,
          dofIndex_to_localIndex),
      _process_data(process_data),
      _e_center_coords(e.getCenterOfGravity().getCoords())
{
    assert(e.getDimension() == GlobalDim - 1);

    IntegrationMethod integration_method(integration_order);
    unsigned const n_integration_points =
        integration_method.getNumberOfPoints();

    _ip_data.resize(n_integration_points);

    auto const shape_matrices_p =
        initShapeMatrices<ShapeFunctionPressure, ShapeMatricesTypePressure,
                          IntegrationMethod, GlobalDim>(e, is_axially_symmetric,
                                                        integration_method);

    auto mat_id = (*_process_data.mesh_prop_materialIDs)[e.getID()];
    auto frac_id = _process_data.map_materialID_to_fractureID[mat_id];
    _fracture_property = &*_process_data.fracture_properties[frac_id];

    // collect properties of connecting fractures, junctions
    for (auto fid : process_data.vec_ele_connected_fractureIDs[e.getID()])
    {
        _fracID_to_local.insert({fid, _fracture_props.size()});
        _fracture_props.push_back(&*_process_data.fracture_properties[fid]);
    }
    for (auto jid : process_data.vec_ele_connected_junctionIDs[e.getID()])
        _junction_props.push_back(&_process_data.junction_properties[jid]);

    typename ShapeMatricesTypePressure::NodalVectorType
        aperture0_node_values =
            _fracture_property->aperture0.getNodalValuesOnElement(
            e, /*time independent*/ 0);

    ParameterLib::SpatialPosition x_position;
    x_position.setElementID(e.getID());
    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        x_position.setIntegrationPoint(ip);

        auto const& sm_p = shape_matrices_p[ip];
        auto& ip_data = _ip_data[ip];
        ip_data.integration_weight =
            sm_p.detJ * sm_p.integralMeasure *
            integration_method.getWeightedPoint(ip).getWeight();

        ip_data.N_p = sm_p.N;
        ip_data.dNdx_p = sm_p.dNdx;

        ip_data.aperture0 = aperture0_node_values.dot(sm_p.N);
        ip_data.aperture = ip_data.aperture0;

        ip_data.permeability_state =
            _fracture_property->permeability_model->getNewState();
    }
}

template <typename ShapeFunctionPressure, typename IntegrationMethod,
          int GlobalDim>
void THLocalAssemblerFracture<
    ShapeFunctionPressure, IntegrationMethod,
    GlobalDim>::assembleWithJacobianConcrete(double const t,
                                             Eigen::VectorXd const& local_x,
                                             Eigen::VectorXd const& local_xdot,
                                             Eigen::VectorXd& local_b,
                                             Eigen::MatrixXd& local_J)
{
    // auto const n_fractures = _fracture_props.size();
    // auto const n_junctions = _junction_props.size();
    // auto const n_enrich_var = n_fractures + n_junctions;

    auto const local_p = local_x.segment<pressure_size>(pressure_index);
    auto const local_p_dot = local_xdot.segment<pressure_size>(pressure_index);
    auto const local_T = local_x.segment<temperature_size>(temperature_index);
    auto const local_T_dot =
        local_xdot.segment<temperature_size>(temperature_index);

    //
    auto local_b_p = local_b.segment<pressure_size>(pressure_index);
    auto local_b_T = local_b.segment<temperature_size>(temperature_index);
    auto local_J_pp = local_J.block<pressure_size, pressure_size>(
        pressure_index, pressure_index);
    auto local_J_pT = local_J.block<pressure_size, temperature_size>(
        pressure_index, temperature_index);
    auto local_J_Tp = local_J.block<temperature_size, pressure_size>(
        temperature_index, pressure_index);
    auto local_J_TT = local_J.block<temperature_size, temperature_size>(
        temperature_index, temperature_index);

    assembleBlockMatricesWithJacobian(
        t, local_p, local_p_dot, local_T, local_T_dot, local_b_p, local_b_T,
        local_J_pp, local_J_pT, local_J_TT, local_J_Tp);
}

template <typename ShapeFunctionPressure, typename IntegrationMethod,
          int GlobalDim>
void THLocalAssemblerFracture<ShapeFunctionPressure, IntegrationMethod,
                              GlobalDim>::
    assembleBlockMatricesWithJacobian(
        double const t, Eigen::Ref<const Eigen::VectorXd> const& p,
        Eigen::Ref<const Eigen::VectorXd> const& p_dot,
        Eigen::Ref<const Eigen::VectorXd> const& T,
        Eigen::Ref<const Eigen::VectorXd> const& T_dot,
        Eigen::Ref<Eigen::VectorXd> rhs_p, Eigen::Ref<Eigen::VectorXd> rhs_T,
        Eigen::Ref<Eigen::MatrixXd> J_pp, Eigen::Ref<Eigen::MatrixXd> /*J_pT*/,
        Eigen::Ref<Eigen::MatrixXd> J_TT, Eigen::Ref<Eigen::MatrixXd> /*J_Tp*/
    )
{
    auto const& R = _fracture_property->R;
    double const& dt = _process_data.dt;

    using GlobalDimMatrix = Eigen::Matrix<double, GlobalDim, GlobalDim>;
    // using GlobalDimVector = Eigen::Matrix<double, GlobalDim, 1>;

    auto const& b = _process_data.specific_body_force;

    ParameterLib::SpatialPosition x_position;
    x_position.setElementID(_element.getID());

    unsigned const n_integration_points = _ip_data.size();
    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        x_position.setIntegrationPoint(ip);

        auto& ip_data = _ip_data[ip];
        auto const& ip_w = ip_data.integration_weight;
        auto const& N_p = ip_data.N_p;
        auto const& dNdx_p = ip_data.dNdx_p;
        auto const& N_T = ip_data.N_p;
        auto const& dNdx_T = ip_data.dNdx_p;

        // auto& mat = ip_data.fracture_material;

        auto const p_dot_ip = N_p * p_dot;
        // auto const dp_ip = p_dot_ip * dt;
        double const p1_ip = N_p * p;
        // auto const p0_ip = p1_ip - dp_ip;
        auto const grad_p1 = (dNdx_p * p).eval();

        auto const T_dot_ip = N_T * T_dot;
        // auto const dT_ip = T_dot_ip * dt;
        auto const T1_ip = N_T * T;
        // auto const T0_ip = T1_ip - dT_ip;
        auto const grad_T1 = (dNdx_T * T).eval();

        auto& b_m = ip_data.aperture0;
        auto& q = ip_data.q;

        //------------------------------------------------------
        // fluid properties
        //------------------------------------------------------
        using ArrayType = MaterialLib::Fluid::FluidProperty::ArrayType;
        using VariableType = MaterialLib::Fluid::PropertyVariableType;
        ArrayType fluid_vars;
        fluid_vars[static_cast<int>(VariableType::T)] = T1_ip;
        fluid_vars[static_cast<int>(VariableType::p)] = p1_ip;

        auto const rho_f = _process_data.fluid_density.getValue(fluid_vars);
        auto const drhof_dp =
            _process_data.fluid_density.getdValue(fluid_vars, VariableType::p);
        auto const drhof_dT =
            _process_data.fluid_density.getdValue(fluid_vars, VariableType::T);
        auto const beta_T_f = drhof_dT / rho_f;  // TODO
        double const mu = _process_data.fluid_viscosity.getValue(fluid_vars);
        auto const dmu_dp = _process_data.fluid_viscosity.getdValue(
            fluid_vars, VariableType::p);
        auto const dmu_dT = _process_data.fluid_viscosity.getdValue(
            fluid_vars, VariableType::T);
        double const lambda_f =
            _process_data.fluid_thermal_conductivity.getValue(fluid_vars);
        double const cp_f =
            _process_data.fluid_specific_heat_capacity.getValue(fluid_vars);
        double const dcpf_dp =
            _process_data.fluid_specific_heat_capacity.getdValue(
                fluid_vars, VariableType::p);
        double const dcpf_dT =
            _process_data.fluid_specific_heat_capacity.getdValue(
                fluid_vars, VariableType::T);
        double const Cp_f = rho_f * cp_f;
        ;
        double const dCpf_dp = drhof_dp * cp_f + rho_f * dcpf_dp;
        double const dCpf_dT = drhof_dT * cp_f + rho_f * dcpf_dT;

        //------------------------------------------------------
        // bulk properties
        //------------------------------------------------------
        double const S =
            (*_fracture_property->specific_storage)(t, x_position)[0];
        // double const biot =
        //     (*_fracture_property->biot_coefficient)(t, x_position)[0];
        auto& permeability = ip_data.permeability;
        permeability = _fracture_property->permeability_model->permeability(
            ip_data.permeability_state.get(), ip_data.aperture0, b_m);

        GlobalDimMatrix const k =
            createRotatedTensor<GlobalDim>(R, permeability);

        //------------------------------------------------------
        // Darcy flux calculation
        //------------------------------------------------------
        q.noalias() = -k / mu * (grad_p1 - rho_f * b);
        auto dq_dpi = (-k / mu * dNdx_p).eval();
        dq_dpi.noalias() += k / mu * drhof_dp * b * N_p;
        dq_dpi.noalias() += k / mu / mu * dmu_dp * (grad_p1 - rho_f * b) * N_p;
        auto dq_dTi = (k / mu * drhof_dT * b * N_T).eval();
        dq_dTi.noalias() += k / mu / mu * dmu_dT * (grad_p1 - rho_f * b) * N_T;

        //------------------------------------------------------
        // Heat flux calculation
        //------------------------------------------------------
        // conduction
        auto const jdiff = (-lambda_f * grad_T1).eval();
        auto const djdiff_dTi = (-lambda_f * dNdx_T).eval();
        // advection
        auto const djadv_dx = (q.transpose() * Cp_f * grad_T1).eval();
        auto ddjadvdx_dpi = (q.transpose() * dCpf_dp * grad_T1 * N_p).eval();
        ddjadvdx_dpi.noalias() += dq_dpi.transpose() * Cp_f * grad_T1;
        auto ddjadvdx_dTi = (q.transpose() * Cp_f * dNdx_T).eval();
        ddjadvdx_dTi.noalias() += q.transpose() * dCpf_dT * grad_T1 * N_T;
        ddjadvdx_dTi.noalias() += dq_dTi.transpose() * Cp_f * grad_T1;

        //------------------------------------------------------
        // residual calculations
        //------------------------------------------------------
        rhs_p.noalias() -= N_p.transpose() * b_m *
                               (S * p_dot_ip - beta_T_f * T_dot_ip) * ip_w +
                           -dNdx_p.transpose() * b_m * q * ip_w;

        rhs_T.noalias() -= b_m * N_T.transpose() * Cp_f * T_dot_ip * ip_w -
                           b_m * dNdx_T.transpose() * jdiff * ip_w +
                           b_m * N_T.transpose() * djadv_dx * ip_w;

        //------------------------------------------------------
        // jacobian calculations
        //------------------------------------------------------
        J_pp.noalias() += b_m * N_p.transpose() * S * N_p * ip_w / dt;
        J_pp.noalias() += -b_m * dNdx_p.transpose() * dq_dpi * ip_w;
        // J_pT.noalias() += - b_m * N_p.transpose() * beta_T_f * N_T * ip_w /
        // dt; J_pT.noalias() += - b_m * dNdx_p.transpose() * dq_dTi * ip_w;

        J_TT.noalias() += b_m * N_T.transpose() * Cp_f * N_T * ip_w / dt;
        // J_TT.noalias() += b_m * N_T.transpose() * dCpf_dT * T_dot_ip * N_T *
        // ip_w;
        J_TT.noalias() += -b_m * dNdx_T.transpose() * djdiff_dTi * ip_w;
        // J_Tp.noalias() += b_m * N_T.transpose() * dCpf_dp * T_dot_ip * N_p *
        // ip_w;
        J_TT.noalias() += b_m * N_T.transpose() * ddjadvdx_dTi * ip_w;
    }

    // std::cout << "Element: " << _element.getID() << "\n";
    // std::cout << "p:\n" << p << "\n";
    // std::cout << "g:\n" << g << "\n";
    // std::cout << "r_p:\n" << rhs_p << "\n";
    // std::cout << "r_g:\n" << rhs_g << "\n";
    // std::cout << "J_pp:\n" << J_pp << "\n";
    // std::cout << "J_pg:\n" << J_pg << "\n";
    // std::cout << "J_gp:\n" << J_gp << "\n";
    // std::cout << "J_gu:\n" << J_gg << "\n";
}

template <typename ShapeFunctionPressure, typename IntegrationMethod,
          int GlobalDim>
void THLocalAssemblerFracture<ShapeFunctionPressure, IntegrationMethod,
                              GlobalDim>::
    computeSecondaryVariableConcreteWithVector(const double /*t*/,
                                               Eigen::VectorXd const& /*local_x*/)
{
    // auto const n_fractures = _fracture_props.size();
    // auto const n_junctions = _junction_props.size();
    // auto const n_enrich_var = n_fractures + n_junctions;

    // auto const& R = _fracture_property->R;

    ParameterLib::SpatialPosition x_position;
    x_position.setElementID(_element.getID());

    unsigned const n_integration_points = _ip_data.size();

    Eigen::Vector3d ele_velocity = Eigen::Vector3d::Zero();
    for (auto const& ip : _ip_data)
    {
        ele_velocity.head(GlobalDim) += ip.q;
    }
    ele_velocity /= static_cast<double>(n_integration_points);
    auto const element_id = _element.getID();
    for (unsigned i = 0; i < 3; i++)
    {
        (*_process_data.mesh_prop_velocity)[element_id * 3 + i] =
            ele_velocity[i];
    }
}

}  // namespace TH
}  // namespace LIE
}  // namespace ProcessLib
