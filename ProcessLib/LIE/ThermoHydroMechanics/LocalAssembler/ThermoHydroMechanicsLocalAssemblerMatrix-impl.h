
#pragma once

#include <iostream>

#include "ThermoHydroMechanicsLocalAssemblerMatrix.h"

#include "MaterialLib/PhysicalConstant.h"
#include "MaterialLib/SolidModels/SelectSolidConstitutiveRelation.h"
#include "MathLib/KelvinVector.h"
#include "MeshLib/ElementStatus.h"
#include "NumLib/Function/Interpolation.h"
#include "ProcessLib/Deformation/LinearBMatrix.h"
#include "ProcessLib/Utils/InitShapeMatrices.h"

namespace ProcessLib
{
namespace LIE
{
namespace ThermoHydroMechanics
{
template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int GlobalDim>
ThermoHydroMechanicsLocalAssemblerMatrix<ShapeFunctionDisplacement,
                                   ShapeFunctionPressure, IntegrationMethod,
                                   GlobalDim>::
    ThermoHydroMechanicsLocalAssemblerMatrix(
        MeshLib::Element const& e,
        std::size_t const n_variables,
        std::size_t const /*local_matrix_size*/,
        std::vector<unsigned> const& dofIndex_to_localIndex,
        bool const is_axially_symmetric,
        unsigned const integration_order,
        ThermoHydroMechanicsProcessData<GlobalDim>& process_data)
    : ThermoHydroMechanicsLocalAssemblerInterface(
          e, is_axially_symmetric,
          (n_variables - 2) * ShapeFunctionDisplacement::NPOINTS * GlobalDim +
              ShapeFunctionPressure::NPOINTS * 2,
          dofIndex_to_localIndex),
      _process_data(process_data)
{
    IntegrationMethod integration_method(integration_order);
    unsigned const n_integration_points =
        integration_method.getNumberOfPoints();

    _ip_data.reserve(n_integration_points);

    auto const shape_matrices_u =
        initShapeMatrices<ShapeFunctionDisplacement,
                          ShapeMatricesTypeDisplacement, IntegrationMethod,
                          GlobalDim>(e, is_axially_symmetric,
                                     integration_method);

    auto const shape_matrices_p =
        initShapeMatrices<ShapeFunctionPressure, ShapeMatricesTypePressure,
                          IntegrationMethod, GlobalDim>(e, is_axially_symmetric,
                                                        integration_method);

    auto& solid_material = MaterialLib::Solids::selectSolidConstitutiveRelation(
        _process_data.solid_materials, _process_data.material_ids, e.getID());

    ParameterLib::SpatialPosition x_position;
    x_position.setElementID(e.getID());
    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        x_position.setIntegrationPoint(ip);

        _ip_data.emplace_back(solid_material);
        auto& ip_data = _ip_data[ip];
        auto const& sm_u = shape_matrices_u[ip];
        auto const& sm_p = shape_matrices_p[ip];
        ip_data.integration_weight =
            sm_u.detJ * sm_u.integralMeasure *
            integration_method.getWeightedPoint(ip).getWeight();

        ip_data.N_u = sm_u.N;
        ip_data.dNdx_u = sm_u.dNdx;
        ip_data.H_u.setZero(GlobalDim, displacement_size);
        for (int i = 0; i < GlobalDim; ++i)
        {
            ip_data.H_u
                .template block<1, displacement_size / GlobalDim>(
                    i, i * displacement_size / GlobalDim)
                .noalias() = ip_data.N_u;
        }

        ip_data.N_p = sm_p.N;
        ip_data.dNdx_p = sm_p.dNdx;

        ip_data.sigma_eff.setZero(kelvin_vector_size);
        ip_data.eps.setZero(kelvin_vector_size);

        ip_data.sigma_eff_prev.resize(kelvin_vector_size);
        ip_data.eps_prev.resize(kelvin_vector_size);
        ip_data.C.resize(kelvin_vector_size, kelvin_vector_size);

        auto const initial_effective_stress =
            _process_data.initial_effective_stress(0, x_position);
        for (unsigned i = 0; i < kelvin_vector_size; i++)
        {
            ip_data.sigma_eff[i] = initial_effective_stress[i];
            ip_data.sigma_eff_prev[i] = initial_effective_stress[i];
        }
    }
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int GlobalDim>
void ThermoHydroMechanicsLocalAssemblerMatrix<
    ShapeFunctionDisplacement, ShapeFunctionPressure, IntegrationMethod,
    GlobalDim>::assembleWithJacobianConcrete(double const t,
                                             Eigen::VectorXd const& local_x,
                                             Eigen::VectorXd const& local_x_dot,
                                             Eigen::VectorXd& local_rhs,
                                             Eigen::MatrixXd& local_Jac)
{
    auto p = const_cast<Eigen::VectorXd&>(local_x).segment(pressure_index,
                                                           pressure_size);
    auto p_dot = const_cast<Eigen::VectorXd&>(local_x_dot)
                     .segment(pressure_index, pressure_size);

    if (_process_data.deactivate_matrix_in_flow)
    {
        setPressureOfInactiveNodes(t, p);
        setPressureDotOfInactiveNodes(p_dot);
    }

    auto T = const_cast<Eigen::VectorXd&>(local_x).segment(temperature_index,
                                                           temperature_size);
    auto T_dot = const_cast<Eigen::VectorXd&>(local_x_dot)
                     .segment(temperature_index, temperature_size);

    auto u = local_x.segment(displacement_index, displacement_size);
    auto u_dot = local_x_dot.segment(displacement_index, displacement_size);

    auto rhs_p = local_rhs.template segment<pressure_size>(pressure_index);
    auto rhs_T = local_rhs.template segment<temperature_size>(temperature_index);
    auto rhs_u =
        local_rhs.template segment<displacement_size>(displacement_index);

    auto J_pp = local_Jac.template block<pressure_size, pressure_size>(
        pressure_index, pressure_index);
    auto J_pT = local_Jac.template block<pressure_size, temperature_size>(
        pressure_index, temperature_index);
    auto J_pu = local_Jac.template block<pressure_size, displacement_size>(
        pressure_index, displacement_index);
    auto J_TT = local_Jac.template block<temperature_size, temperature_size>(
        temperature_index, temperature_index);
    auto J_Tp = local_Jac.template block<temperature_size, pressure_size>(
        temperature_index, pressure_index);
    auto J_Tu = local_Jac.template block<temperature_size, displacement_size>(
        temperature_index, displacement_index);
    auto J_uu = local_Jac.template block<displacement_size, displacement_size>(
        displacement_index, displacement_index);
    auto J_up = local_Jac.template block<displacement_size, pressure_size>(
        displacement_index, pressure_index);
    auto J_uT = local_Jac.template block<displacement_size, temperature_size>(
        displacement_index, temperature_index);

    assembleBlockMatricesWithJacobian(t, p, p_dot, T, T_dot, u, u_dot, rhs_p, rhs_T, rhs_u, J_pp,
                                      J_pT, J_pu, J_TT, J_Tp, J_Tu, J_uu, J_up, J_uT);
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int GlobalDim>
void ThermoHydroMechanicsLocalAssemblerMatrix<ShapeFunctionDisplacement,
                                              ShapeFunctionPressure,
                                              IntegrationMethod, GlobalDim>::
    assembleBlockMatricesWithJacobian(
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
            /*J_pT*/,
        Eigen::Ref<Eigen::MatrixXd>
            J_pu,
        Eigen::Ref<Eigen::MatrixXd>
            J_TT,
        Eigen::Ref<Eigen::MatrixXd>
            /*J_Tp*/,
        Eigen::Ref<Eigen::MatrixXd>
            /*J_Tu*/,
        Eigen::Ref<Eigen::MatrixXd>
            J_uu,
        Eigen::Ref<Eigen::MatrixXd>
            J_up,
        Eigen::Ref<Eigen::MatrixXd>
            /*J_uT*/)
{
    assert(this->_element.getDimension() == GlobalDim);

    double const& dt = _process_data.dt;
    auto const& b = _process_data.specific_body_force;

    ParameterLib::SpatialPosition x_position;
    x_position.setElementID(_element.getID());

    using Invariants = MathLib::KelvinVector::Invariants<
        MathLib::KelvinVector::KelvinVectorDimensions<GlobalDim>::value>;

    unsigned const n_integration_points = _ip_data.size();
    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        x_position.setIntegrationPoint(ip);

        auto& ip_data = _ip_data[ip];
        auto const& ip_w = ip_data.integration_weight;
        auto const& N_u = ip_data.N_u;
        auto const& dNdx_u = ip_data.dNdx_u;
        auto const& N_p = ip_data.N_p;
        auto const& dNdx_p = ip_data.dNdx_p;
        auto const& H_u = ip_data.H_u;
        auto const& N_T = N_p;
        auto const& dNdx_T = dNdx_p;

        auto const x_coord =
            interpolateXCoordinate<ShapeFunctionDisplacement,
                                   ShapeMatricesTypeDisplacement>(_element,
                                                                  N_u);
        auto const B =
            LinearBMatrix::computeBMatrix<GlobalDim,
                                          ShapeFunctionDisplacement::NPOINTS,
                                          typename BMatricesType::BMatrixType>(
                dNdx_u, N_u, x_coord, _is_axially_symmetric);

        auto const p_dot_ip = N_p.dot(p_dot);
        // auto const dp_ip = p_dot_ip * dt;
        auto const p1_ip = N_p * p;
        // auto const p0_ip = p1_ip - dp_ip;
        auto const grad_p1_ip = (dNdx_p * p).eval();
        auto const T_dot_ip = N_T.dot(T_dot);
        auto const dT_ip = T_dot_ip * dt;
        auto const T1_ip = N_T * T;
        //auto const T0_ip = T1_ip - dT_ip;
        auto const grad_T1_ip = (dNdx_T * T).eval();

        auto& eps = ip_data.eps;
        auto const& eps_prev = ip_data.eps_prev;
        auto& eps_m = ip_data.eps_m;
        auto const& eps_m_prev = ip_data.eps_m_prev;
        auto& sigma = ip_data.sigma;
        // auto const& sigma_prev = ip_data.sigma_prev;
        auto& sigma_eff = ip_data.sigma_eff;
        auto& sigma_eff_prev = ip_data.sigma_eff_prev;
        auto& q = ip_data.q;

        auto& state = ip_data.material_state_variables;

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

        //------------------------------------------------------
        // solid properties
        //------------------------------------------------------
        auto const rho_s = _process_data.solid_density(t, x_position)[0];
        // auto const drhos_dp = 0.0;
        // auto const drhos_dT = 0.0;
        //double const rho_s = rho_sr * (1 - 3 * thermal_strain);
        double const alpha_T_s =
            _process_data.solid_linear_thermal_expansion_coefficient(
                t, x_position)[0];
        double const beta_T_s = 3 * alpha_T_s;
        double const lambda_s =
            _process_data.solid_thermal_conductivity(t, x_position)[0];
        double const cp_s =
            _process_data.solid_specific_heat_capacity(t, x_position)[0];

        //------------------------------------------------------
        // bulk properties
        //------------------------------------------------------
        double const S = _process_data.specific_storage(t, x_position)[0];
        auto const ki = _process_data.intrinsic_permeability(t, x_position)[0];
        auto const biot = _process_data.biot_coefficient(t, x_position)[0];
        auto const porosity = _process_data.porosity(t, x_position)[0];
        auto const rho = rho_s * (1. - porosity) + porosity * rho_f;
        // auto const drho_dp = drhos_dp * (1. - porosity) + porosity * drhof_dp;
        // auto const drho_dT = drhos_dT * (1. - porosity) + porosity * drhof_dT;
        auto const lambda = porosity * lambda_f + (1 - porosity) * lambda_s;
/* Lauwerier benchmark needs anisotropic tensor
        Eigen::MatrixXd lambda(GlobalDim, GlobalDim);
        lambda.setZero();
        lambda(1,1) = (1 - porosity) * lambda_s;
*/
        auto const Cp =
            porosity * cp_f * rho_f + (1 - porosity) * cp_s * rho_s;
        // auto const dCp_dp =
        //     porosity * cp_f * drhof_dp + (1 - porosity) * cp_s * drhos_dp;
        // auto const dCp_dT =
        //     porosity * cp_f * drhof_dT + (1 - porosity) * cp_s * drhos_dT;
        auto const beta_T = porosity * beta_T_f + (1 - porosity) * beta_T_s;

        //------------------------------------------------------
        // Darcy flux calculation
        //------------------------------------------------------
        q.noalias() = -ki / mu * (grad_p1_ip - rho_f * b);
        auto dq_dpi = (-ki / mu * dNdx_p).eval();
        dq_dpi.noalias() += ki / mu * drhof_dp * b * N_p;
        dq_dpi.noalias() += ki / mu / mu * dmu_dp * grad_p1_ip * N_p;
        auto dq_dTi = (ki / mu * drhof_dT * b * N_T).eval();
        dq_dTi.noalias() += ki / mu / mu * dmu_dT * grad_p1_ip * N_T;
        if (_process_data.deactivate_matrix_in_flow)
        {
            q.setZero();
            dq_dpi.setZero();
            dq_dTi.setZero();
        }

        //------------------------------------------------------
        // Heat flux calculation
        //------------------------------------------------------
        // conduction
        auto const jdiff = (- lambda * grad_T1_ip).eval();
        auto const djdiff_dTi = (- lambda * dNdx_T).eval();
        // advection
        auto djadv_dx = (q.transpose() * rho_f * cp_f * grad_T1_ip).eval();
        auto ddjadvdx_dpi = (q.transpose() * (drhof_dp * cp_f + rho_f * dcpf_dp ) * grad_T1_ip * N_p).eval();
        ddjadvdx_dpi.noalias() += dq_dpi.transpose() * rho_f * cp_f * grad_T1_ip;
        auto ddjadvdx_dTi = (q.transpose() * rho_f * cp_f * dNdx_T).eval();
        ddjadvdx_dTi.noalias() +=
            q.transpose() * (drhof_dT * cp_f + rho_f * dcpf_dT ) * grad_T1_ip * N_T;
        ddjadvdx_dTi.noalias() += dq_dTi.transpose() * rho_f * cp_f * grad_T1_ip;
        if (_process_data.deactivate_matrix_in_flow)
        {
            djadv_dx.setZero();
            ddjadvdx_dpi.setZero();
            ddjadvdx_dTi.setZero();
        }

        //------------------------------------------------------
        // strain calculation
        //------------------------------------------------------
        eps.noalias() = B * u;
        auto const eps_T = alpha_T_s * dT_ip;
        eps_m.noalias() =
            eps_m_prev + eps - eps_prev - eps_T * Invariants::identity2;
        auto const eps_dot = B * u_dot;

        //------------------------------------------------------
        // stress, C calculation
        //------------------------------------------------------
        // if (sigma_prev.isZero(0)) {
        //     sigma_eff_prev.noalias() = sigma_prev;
        // } else {
        //     sigma_eff_prev.noalias() = sigma_prev + biot * p0_ip * Invariants::identity2;
        // }
        auto&& solution = _ip_data[ip].solid_material.integrateStress(
            t, x_position, dt, eps_m_prev, eps_m, sigma_eff_prev, *state, T1_ip);

        if (!solution)
            OGS_FATAL("Computation of local constitutive relation failed.");

        MathLib::KelvinVector::KelvinMatrixType<GlobalDim> C;
        std::tie(sigma_eff, state, C) = std::move(*solution);
        sigma.noalias() = sigma_eff - biot * Invariants::identity2 * p1_ip;

        //------------------------------------------------------
        // residual calculations
        //------------------------------------------------------
        // pressure equation
        if (!_process_data.deactivate_matrix_in_flow)
        {
            rhs_p.noalias() -=
                N_p.transpose() *
                    (S * p_dot_ip - beta_T * T_dot_ip +
                     biot * Invariants::identity2.transpose() * eps_dot) *
                    ip_w
                - dNdx_p.transpose() * q * ip_w;
        }

        // temperature equation
        rhs_T.noalias() -=
            N_T.transpose() * Cp * T_dot_ip * ip_w
            - dNdx_T.transpose() * jdiff * ip_w;
        if (!_process_data.deactivate_matrix_in_flow)
            rhs_T.noalias() -= N_T.transpose() * djadv_dx * ip_w;

        // displacement equation
        rhs_u.noalias() -=
            B.transpose() * sigma * ip_w
            - H_u.transpose() * rho * b * ip_w;


        //------------------------------------------------------
        // jacobian calculations
        //------------------------------------------------------
        if (!_process_data.deactivate_matrix_in_flow)
        {
            J_pp.noalias() += N_p.transpose() * S * N_p * ip_w / dt;
            J_pp.noalias() += -dNdx_p.transpose() * dq_dpi * ip_w;
            // J_pT.noalias() += -N_p.transpose() * beta_T * N_T * ip_w / dt;
            // J_pT.noalias() += -dNdx_p.transpose() * dq_dTi * ip_w;
            J_pu.noalias() += N_p.transpose() * biot *
                              Invariants::identity2.transpose() * B * ip_w / dt;
        }

        J_TT.noalias() += N_T.transpose() * Cp * N_T * ip_w / dt;
        // J_TT.noalias() += N_T.transpose() * dCp_dT * T_dot_ip * N_T * ip_w;
        J_TT.noalias() += - dNdx_T.transpose() * djdiff_dTi * ip_w;
        // J_Tp.noalias() += N_T.transpose() * dCp_dp * T_dot_ip * N_p * ip_w;
        // if (!_process_data.deactivate_matrix_in_flow)
        // {
        //     J_TT.noalias() += N_T.transpose() * ddjadvdx_dTi * ip_w;
        //     J_Tp.noalias() += N_T.transpose() * ddjadvdx_dpi * ip_w;
        // }

        J_uu.noalias() += B.transpose() * C * B * ip_w;
        J_up.noalias() += - B.transpose() * biot * Invariants::identity2 * N_p * ip_w;
        // J_up.noalias() += - H_u.transpose() * drho_dp * b * N_p * ip_w;
        // J_uT.noalias() +=
        //     B.transpose() * C * alpha_T_s * Invariants::identity2 * N_T * ip_w;
        // J_uT.noalias() += - H_u.transpose() * drho_dT * b * N_T * ip_w;
    }

    // std::cout << "Element: " << _element.getID() << "\n";
    // std::cout << "r_p:\n" << rhs_p << "\n"; 
    // std::cout << "r_u:\n" << rhs_u << "\n"; 
    // std::cout << "J_pp:\n" << J_pp << "\n"; 
    // std::cout << "J_pu:\n" << J_pu << "\n"; 
    // std::cout << "J_up:\n" << J_up << "\n"; 
    // std::cout << "J_uu:\n" << J_uu << "\n"; 
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int GlobalDim>
void ThermoHydroMechanicsLocalAssemblerMatrix<ShapeFunctionDisplacement,
                                        ShapeFunctionPressure,
                                        IntegrationMethod, GlobalDim>::
    computeSecondaryVariableConcreteWithVector(double const t,
                                               Eigen::VectorXd const& local_x)
{
    auto p = const_cast<Eigen::VectorXd&>(local_x).segment(pressure_index,
                                                           pressure_size);
    if (_process_data.deactivate_matrix_in_flow)
    {
        setPressureOfInactiveNodes(t, p);
    }
    auto T = const_cast<Eigen::VectorXd&>(local_x).segment(temperature_index,
        temperature_size);
    auto u = local_x.segment(displacement_index, displacement_size);

    computeSecondaryVariableConcreteWithBlockVectors(t, p, T, u);
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int GlobalDim>
void ThermoHydroMechanicsLocalAssemblerMatrix<ShapeFunctionDisplacement,
                                        ShapeFunctionPressure,
                                        IntegrationMethod, GlobalDim>::
    computeSecondaryVariableConcreteWithBlockVectors(
        double const /*t*/,
        Eigen::Ref<const Eigen::VectorXd> const& p,
        Eigen::Ref<const Eigen::VectorXd> const& T,
        Eigen::Ref<const Eigen::VectorXd> const& /*u*/)
{
    int n = GlobalDim == 2 ? 4 : 6;
    Eigen::VectorXd ele_stress = Eigen::VectorXd::Zero(n);
    Eigen::VectorXd ele_strain = Eigen::VectorXd::Zero(n);
    Eigen::Vector3d ele_velocity = Eigen::Vector3d::Zero();

    auto const n_integration_points = _ip_data.size();
    for (auto const& ip_data : _ip_data)
    {
        ele_stress += ip_data.sigma_eff;
        ele_strain += ip_data.eps;
        ele_velocity.head(GlobalDim) += ip_data.q;
    }

    ele_stress /= static_cast<double>(n_integration_points);
    ele_strain /= static_cast<double>(n_integration_points);
    ele_velocity /= static_cast<double>(n_integration_points);

    auto const element_id = _element.getID();
    (*_process_data.mesh_prop_stress_xx)[_element.getID()] = ele_stress[0];
    (*_process_data.mesh_prop_stress_yy)[_element.getID()] = ele_stress[1];
    (*_process_data.mesh_prop_stress_zz)[_element.getID()] = ele_stress[2];
    (*_process_data.mesh_prop_stress_xy)[_element.getID()] = ele_stress[3];
    if (GlobalDim == 3)
    {
        (*_process_data.mesh_prop_stress_yz)[_element.getID()] = ele_stress[4];
        (*_process_data.mesh_prop_stress_xz)[_element.getID()] = ele_stress[5];
    }

    (*_process_data.mesh_prop_strain_xx)[_element.getID()] = ele_strain[0];
    (*_process_data.mesh_prop_strain_yy)[_element.getID()] = ele_strain[1];
    (*_process_data.mesh_prop_strain_zz)[_element.getID()] = ele_strain[2];
    (*_process_data.mesh_prop_strain_xy)[_element.getID()] = ele_strain[3];
    if (GlobalDim == 3)
    {
        (*_process_data.mesh_prop_strain_yz)[_element.getID()] = ele_strain[4];
        (*_process_data.mesh_prop_strain_xz)[_element.getID()] = ele_strain[5];
    }

    for (unsigned i = 0; i < 3; i++)
    {
        (*_process_data.mesh_prop_velocity)[element_id * 3 + i] =
            ele_velocity[i];
    }

    NumLib::interpolateToHigherOrderNodes<
        ShapeFunctionPressure, typename ShapeFunctionDisplacement::MeshElement,
        GlobalDim>(_element, _is_axially_symmetric, p,
                   *_process_data.mesh_prop_nodal_p);

    NumLib::interpolateToHigherOrderNodes<
        ShapeFunctionPressure, typename ShapeFunctionDisplacement::MeshElement,
        GlobalDim>(_element, _is_axially_symmetric, T,
                   *_process_data.mesh_prop_nodal_T);
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int GlobalDim>
void ThermoHydroMechanicsLocalAssemblerMatrix<
    ShapeFunctionDisplacement, ShapeFunctionPressure, IntegrationMethod,
    GlobalDim>::setPressureOfInactiveNodes(double const t,
                                           Eigen::Ref<Eigen::VectorXd>
                                               p)
{
    ParameterLib::SpatialPosition x_position;
    x_position.setElementID(_element.getID());
    for (unsigned i = 0; i < pressure_size; i++)
    {
        // only inactive nodes
        if (_process_data.p_element_status->isActiveNode(_element.getNode(i)))
        {
            continue;
        }
        x_position.setNodeID(_element.getNodeIndex(i));
        auto const p0 = (*_process_data.p0)(t, x_position)[0];
        p[i] = p0;
    }
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int GlobalDim>
void ThermoHydroMechanicsLocalAssemblerMatrix<
    ShapeFunctionDisplacement, ShapeFunctionPressure, IntegrationMethod,
    GlobalDim>::setPressureDotOfInactiveNodes(Eigen::Ref<Eigen::VectorXd> p_dot)
{
    for (unsigned i = 0; i < pressure_size; i++)
    {
        // only inactive nodes
        if (_process_data.p_element_status->isActiveNode(_element.getNode(i)))
        {
            continue;
        }
        p_dot[i] = 0;
    }
}

}  // namespace ThermoHydroMechanics
}  // namespace LIE
}  // namespace ProcessLib
