/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "THMFEM.h"

#include <iostream>

#include "MathLib/KelvinVector.h"
#include "MaterialLib/Fluid/FluidPropertyHeaders.h"
#include "MaterialLib/SolidModels/SelectSolidConstitutiveRelation.h"
#include "NumLib/Function/Interpolation.h"
#include "ProcessLib/CoupledSolutionsForStaggeredScheme.h"

namespace ProcessLib
{
namespace THM
{
template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int DisplacementDim>
THMLocalAssembler<ShapeFunctionDisplacement,
                                   ShapeFunctionPressure, IntegrationMethod,
                                   DisplacementDim>::
    THMLocalAssembler(
        MeshLib::Element const& e,
        std::size_t const /*local_matrix_size*/,
        bool const is_axially_symmetric,
        unsigned const integration_order,
        THMProcessData<DisplacementDim>& process_data)
    : _process_data(process_data),
      _integration_method(integration_order),
      _element(e),
      _is_axially_symmetric(is_axially_symmetric)
{
    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    _ip_data.reserve(n_integration_points);
    _secondary_data.N_u.resize(n_integration_points);

    auto const shape_matrices_u =
        initShapeMatrices<ShapeFunctionDisplacement,
                          ShapeMatricesTypeDisplacement, IntegrationMethod,
                          DisplacementDim>(e, is_axially_symmetric,
                                           _integration_method);

    auto const shape_matrices_p =
        initShapeMatrices<ShapeFunctionPressure, ShapeMatricesTypePressure,
                          IntegrationMethod, DisplacementDim>(
            e, is_axially_symmetric, _integration_method);

    auto const& solid_material =
        MaterialLib::Solids::selectSolidConstitutiveRelation(
            _process_data.solid_materials, _process_data.material_ids,
            e.getID());

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        _ip_data.emplace_back(solid_material);
        auto& ip_data = _ip_data[ip];
        _ip_data[ip].integration_weight =
            _integration_method.getWeightedPoint(ip).getWeight() *
            shape_matrices_u[ip].integralMeasure * shape_matrices_u[ip].detJ;

        ip_data.N_u_op = ShapeMatricesTypeDisplacement::template MatrixType<
            DisplacementDim, displacement_size>::Zero(DisplacementDim,
                                                      displacement_size);
        for (int i = 0; i < DisplacementDim; ++i)
            ip_data.N_u_op
                .template block<1, displacement_size / DisplacementDim>(
                    i, i * displacement_size / DisplacementDim)
                .noalias() = shape_matrices_u[ip].N;

        ip_data.N_u = shape_matrices_u[ip].N;
        ip_data.dNdx_u = shape_matrices_u[ip].dNdx;

        ip_data.N_p = shape_matrices_p[ip].N;
        ip_data.dNdx_p = shape_matrices_p[ip].dNdx;

        _secondary_data.N_u[ip] = shape_matrices_u[ip].N;
    }
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int DisplacementDim>
std::size_t THMLocalAssembler<
    ShapeFunctionDisplacement, ShapeFunctionPressure, IntegrationMethod,
    DisplacementDim>::setIPDataInitialConditions(std::string const& name,
                                                 double const* values,
                                                 int const integration_order)
{
    if (integration_order !=
        static_cast<int>(_integration_method.getIntegrationOrder()))
    {
        OGS_FATAL(
            "Setting integration point initial conditions; The integration "
            "order of the local assembler for element %d is different from "
            "the integration order in the initial condition.",
            _element.getID());
    }

    if (name == "sigma_ip")
    {
        return setSigma(values);
    }
    else if (name == "epsilon_ip")
    {
        return setEpsilon(values);
    }
    else if (name == "epsilon_m_ip")
    {
        return setEpsilonMechanical(values);
    }

    return 0;
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int DisplacementDim>
void THMLocalAssembler<ShapeFunctionDisplacement,
                                        ShapeFunctionPressure,
                                        IntegrationMethod, DisplacementDim>::
    assembleWithJacobian(double const t, std::vector<double> const& local_x,
                         std::vector<double> const& local_xdot,
                         const double /*dxdot_dx*/, const double /*dx_dx*/,
                         std::vector<double>& /*local_M_data*/,
                         std::vector<double>& /*local_K_data*/,
                         std::vector<double>& local_rhs_data,
                         std::vector<double>& local_Jac_data)
{
    assert(local_x.size() ==
           pressure_size + displacement_size + temperature_size);

    auto T = Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
        temperature_size> const>(local_x.data() + temperature_index,
                                 temperature_size);

    auto p = Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
        pressure_size> const>(local_x.data() + pressure_index, pressure_size);

    auto u =
        Eigen::Map<typename ShapeMatricesTypeDisplacement::template VectorType<
            displacement_size> const>(local_x.data() + displacement_index,
                                      displacement_size);

    auto T_dot =
        Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
            temperature_size> const>(local_xdot.data() + temperature_index,
                                     temperature_size);

    auto p_dot =
        Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
            pressure_size> const>(local_xdot.data() + pressure_index,
                                  pressure_size);
    auto u_dot =
        Eigen::Map<typename ShapeMatricesTypeDisplacement::template VectorType<
            displacement_size> const>(local_xdot.data() + displacement_index,
                                      displacement_size);

    auto local_Jac = MathLib::createZeroedMatrix<
        typename ShapeMatricesTypeDisplacement::template MatrixType<
            temperature_size + displacement_size + pressure_size,
            temperature_size + displacement_size + pressure_size>>(
        local_Jac_data, displacement_size + pressure_size + temperature_size,
        displacement_size + pressure_size + temperature_size);

    auto local_rhs = MathLib::createZeroedVector<
        typename ShapeMatricesTypeDisplacement::template VectorType<
            displacement_size + pressure_size + temperature_size>>(
        local_rhs_data, displacement_size + pressure_size + temperature_size);

    auto rhs_T =
        local_rhs.template segment<temperature_size>(temperature_index);
    auto rhs_p = local_rhs.template segment<pressure_size>(pressure_index);
    auto rhs_u =
        local_rhs.template segment<displacement_size>(displacement_index);

    auto Jac_uu = local_Jac.template block<displacement_size, displacement_size>(
                     displacement_index, displacement_index);
    auto Jac_up = local_Jac.template block<displacement_size, pressure_size>(
                     displacement_index, pressure_index);
    auto Jac_uT = local_Jac.template block<displacement_size, temperature_size>(
                     displacement_index, temperature_index);
    auto Jac_pu = local_Jac.template block<pressure_size, displacement_size>(
                     pressure_index, displacement_index);
    auto Jac_pp = local_Jac.template block<pressure_size, pressure_size>(pressure_index,
                                                               pressure_index);
    auto Jac_pT = local_Jac.template block<pressure_size, temperature_size>(
                     pressure_index, temperature_index);
    // auto Jac_Tu = local_Jac.template block<temperature_size, displacement_size>(
    //                  temperature_index, displacement_index);
    auto Jac_Tp = local_Jac.template block<temperature_size, pressure_size>(
                     temperature_index, pressure_index);
    auto Jac_TT = local_Jac.template block<temperature_size, temperature_size>(
                     temperature_index, temperature_index);

    double const& dt = _process_data.dt;

    SpatialPosition x_position;
    x_position.setElementID(_element.getID());

    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();
    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        x_position.setIntegrationPoint(ip);
        auto const& w = _ip_data[ip].integration_weight;

        auto const& N_u_op = _ip_data[ip].N_u_op;
        auto const& N_u = _ip_data[ip].N_u;
        auto const& dNdx_u = _ip_data[ip].dNdx_u;
        auto const& N_p = _ip_data[ip].N_p;
        auto const& dNdx_p = _ip_data[ip].dNdx_p;
        auto const& N_T = N_p;
        auto const& dNdx_T = dNdx_p;

        auto const x_coord =
            interpolateXCoordinate<ShapeFunctionDisplacement,
                                   ShapeMatricesTypeDisplacement>(_element,
                                                                  N_u);
        auto const B =
            LinearBMatrix::computeBMatrix<DisplacementDim,
                                          ShapeFunctionDisplacement::NPOINTS,
                                          typename BMatricesType::BMatrixType>(
                dNdx_u, N_u, x_coord, _is_axially_symmetric);

        auto const dT_ip = N_T.dot(T_dot) * dt;
        auto const T1_ip = N_T * T;
        //auto const T0_ip = T1_ip - dT_ip;
        auto const dp_ip = N_p.dot(p_dot) * dt;
        auto const p1_ip = N_p * p;
        auto const p0_ip = p1_ip - dp_ip;
        auto const grad_p1 = (dNdx_p * p).eval();
        auto const grad_T1 = (dNdx_T * T).eval();

        auto& eps = _ip_data[ip].eps;
        auto const& eps_prev = _ip_data[ip].eps_prev;
        auto& eps_m = _ip_data[ip].eps_m;
        auto const& eps_m_prev = _ip_data[ip].eps_m_prev;
        auto& sigma = _ip_data[ip].sigma;
        auto const& sigma_prev = _ip_data[ip].sigma_prev;
        auto& sigma_eff = _ip_data[ip].sigma_eff;
        auto& sigma_eff_prev = _ip_data[ip].sigma_eff_prev;
        auto& q = _ip_data[ip].q;

        auto& state = _ip_data[ip].material_state_variables;

        //------------------------------------------------------
        // fluid properties
        //------------------------------------------------------
        using ArrayType = MaterialLib::Fluid::FluidProperty::ArrayType;
        ArrayType fluid_vars;
        fluid_vars[static_cast<int>(
            MaterialLib::Fluid::PropertyVariableType::T)] = T1_ip;
        fluid_vars[static_cast<int>(
            MaterialLib::Fluid::PropertyVariableType::p)] = p1_ip;

        auto const rho_f = _process_data.fluid_density.getValue(fluid_vars);
        auto const drhof_dp = _process_data.fluid_density.getdValue(
            fluid_vars, MaterialLib::Fluid::PropertyVariableType::p);
        auto const drhof_dT = _process_data.fluid_density.getdValue(
            fluid_vars, MaterialLib::Fluid::PropertyVariableType::T);
        auto const beta_f = drhof_dT/rho_f; //TODO
        double const mu_f = _process_data.fluid_viscosity.getValue(fluid_vars);
        auto const dmuf_dp = _process_data.fluid_viscosity.getdValue(
            fluid_vars, MaterialLib::Fluid::PropertyVariableType::p);
        auto const dmuf_dT = _process_data.fluid_viscosity.getdValue(
            fluid_vars, MaterialLib::Fluid::PropertyVariableType::T);
        double const lambda_f =
            _process_data.fluid_thermal_conductivity.getValue(fluid_vars);
        double const C_f =
            _process_data.fluid_specific_heat_capacity.getValue(fluid_vars);
        double const dCf_dp =
            _process_data.fluid_specific_heat_capacity.getdValue(
                fluid_vars, MaterialLib::Fluid::PropertyVariableType::p);
        double const dCf_dT =
            _process_data.fluid_specific_heat_capacity.getdValue(
                fluid_vars, MaterialLib::Fluid::PropertyVariableType::T);

        //------------------------------------------------------
        // solid properties
        //------------------------------------------------------
        auto const rho_s = _process_data.solid_density(t, x_position)[0];
        auto const drhos_dp = 0.0;
        auto const drhos_dT = 0.0;
        //double const rho_s = rho_sr * (1 - 3 * thermal_strain);
        double const alpha_s =
            _process_data.solid_linear_thermal_expansion_coefficient(
                t, x_position)[0];
        double const lambda_s =
            _process_data.solid_thermal_conductivity(t, x_position)[0];
        double const C_s =
            _process_data.solid_specific_heat_capacity(t, x_position)[0];

        //------------------------------------------------------
        // bulk properties
        //------------------------------------------------------
        double const S = _process_data.specific_storage(t, x_position)[0];
        auto const ki = _process_data.intrinsic_permeability(t, x_position)[0];
        auto const alpha = _process_data.biot_coefficient(t, x_position)[0];
        auto const porosity = _process_data.porosity(t, x_position)[0];
        auto const rho = rho_s * (1. - porosity) + porosity * rho_f;
        auto const drho_dp = drhos_dp * (1. - porosity) + porosity * drhof_dp;
        auto const drho_dT = drhos_dT * (1. - porosity) + porosity * drhof_dT;
        auto const lambda = porosity * lambda_f + (1 - porosity) * lambda_s;
        auto const heat_capacity =
            porosity * C_f * rho_f + (1 - porosity) * C_s * rho_s;
        auto const dheat_capacity_dp =
            porosity * C_f * drhof_dp + (1 - porosity) * C_s * drhos_dp;
        auto const dheat_capacity_dT =
            porosity * C_f * drhof_dT + (1 - porosity) * C_s * drhos_dT;
        auto const beta = porosity * beta_f + (1 - porosity) * 3 * alpha_s;

        auto const& b = _process_data.specific_body_force;

        //------------------------------------------------------
        // Darcy q calculation
        //------------------------------------------------------
        q.noalias() = -ki / mu_f * (grad_p1 - rho_f * b);
        auto dq_dpi = (-ki / mu_f * dNdx_p).eval();
        dq_dpi.noalias() += ki / mu_f * drhof_dp * b * N_p;
        dq_dpi.noalias() += ki / mu_f / mu_f * dmuf_dp * grad_p1 * N_p;
        auto dq_dTi = (ki / mu_f * drhof_dT * b * N_T).eval();
        dq_dTi.noalias() += ki / mu_f / mu_f * dmuf_dT * grad_p1 * N_T;

        //------------------------------------------------------
        // Heat flux calculation
        //------------------------------------------------------
        auto const jdiff = (- lambda * grad_T1).eval();
        auto const djdiff_dpi = (0.0 * dNdx_p).eval();
        auto const djdiff_dTi = (- lambda * dNdx_T).eval();
        auto const jadv = (q * rho_f * C_f * T1_ip).eval();
        auto djadv_dpi = (q * (drhof_dp * C_f + rho_f * dCf_dp ) * T1_ip * N_p).eval();
        djadv_dpi.noalias() += dq_dpi * rho_f * C_f * T1_ip;
        auto djadv_dTi = (q * rho_f * C_f * N_T).eval();
        djadv_dTi.noalias() +=
            q * (drhof_dT * C_f + rho_f * dCf_dT ) * T1_ip * N_T;
        djadv_dTi.noalias() += dq_dTi * rho_f * C_f * T1_ip;

        //------------------------------------------------------
        // strain calculation
        //------------------------------------------------------
        eps.noalias() = B * u;

        // calculate thermally induced strain
        // assume isotropic thermal expansion
        double const linear_thermal_strain_increment = alpha * dT_ip;

        using Invariants = MathLib::KelvinVector::Invariants<
            MathLib::KelvinVector::KelvinVectorDimensions<
                DisplacementDim>::value>;

        eps_m.noalias() =
            eps_m_prev + eps - eps_prev -
            linear_thermal_strain_increment * Invariants::identity2;

        //------------------------------------------------------
        // stress, C calculation
        //------------------------------------------------------
        sigma_eff_prev.noalias() = sigma_prev + alpha * p0_ip * Invariants::identity2;
        auto&& solution = _ip_data[ip].solid_material.integrateStress(
            t, x_position, dt, eps_m_prev, eps_m, sigma_eff_prev, *state, T1_ip);

        if (!solution)
            OGS_FATAL("Computation of local constitutive relation failed.");

        MathLib::KelvinVector::KelvinMatrixType<DisplacementDim> C;
        std::tie(sigma_eff, state, C) = std::move(*solution);
        sigma.noalias() = sigma_eff - alpha * Invariants::identity2 * p1_ip;

        //------------------------------------------------------
        // residual calculations
        //------------------------------------------------------
        // pressure equation
        rhs_p.noalias() -=
            N_p.transpose() * S * N_p * p_dot * w
            - N_p.transpose() * beta * N_T * T_dot * w
            + N_p.transpose() * alpha * Invariants::identity2.transpose() * (B * u_dot) * w
            - dNdx_p.transpose() * q * w;

        // temperature equation
        rhs_T.noalias() -=
            N_T.transpose() * heat_capacity * N_T * T_dot * w
            - dNdx_T.transpose() * (jdiff + jadv) * w;

        // displacement equation
        rhs_u.noalias() -=
            B.transpose() * sigma * w
            - N_u_op.transpose() * rho * b * w;


        //------------------------------------------------------
        // jacobian calculations
        //------------------------------------------------------
        Jac_pp.noalias() += N_p.transpose() * S * N_p * w / dt;
        Jac_pT.noalias() += - N_p.transpose() * beta * N_T * w / dt;
        Jac_pu.noalias() += (B.transpose() * alpha * Invariants::identity2 * N_p).transpose() * w / dt;
        Jac_pp.noalias() += - dNdx_p.transpose() * dq_dpi * w;
        Jac_pT.noalias() += - dNdx_p.transpose() * dq_dTi * w;

        Jac_TT.noalias() += N_T.transpose() * heat_capacity * N_T * w / dt;
        Jac_Tp.noalias() += N_T.transpose() * dheat_capacity_dp * N_T * T_dot * N_p * w;
        Jac_TT.noalias() += N_T.transpose() * dheat_capacity_dT * N_T * T_dot * N_T * w;
        Jac_TT.noalias() += - dNdx_T.transpose() * (djdiff_dTi + djadv_dTi) * w;
        Jac_Tp.noalias() += - dNdx_T.transpose() * (djdiff_dpi + djadv_dpi) * w;

        Jac_uu.noalias() += B.transpose() * C * B * w;
        Jac_up.noalias() += - B.transpose() * alpha * Invariants::identity2 * N_p * w;
        Jac_up.noalias() += - N_u_op.transpose() * drho_dp * b * N_p * w;
        Jac_uT.noalias() += - N_u_op.transpose() * drho_dT * b * N_T * w;
    }
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int DisplacementDim>
std::vector<double> const& THMLocalAssembler<
    ShapeFunctionDisplacement, ShapeFunctionPressure, IntegrationMethod,
    DisplacementDim>::getIntPtDarcyVelocity(const double /*t*/,
                                            GlobalVector const&
                                                /*current_solution*/,
                                            NumLib::LocalToGlobalIndexMap const&
                                                /*dof_table*/,
                                            std::vector<double>& cache) const
{
    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    cache.clear();
    auto cache_mat = MathLib::createZeroedMatrix<Eigen::Matrix<
        double, DisplacementDim, Eigen::Dynamic, Eigen::RowMajor>>(
        cache, DisplacementDim, n_integration_points);

    for (unsigned ip = 0; ip < n_integration_points; ++ip)
    {
        cache_mat.col(ip).noalias() = _ip_data[ip].q.head(DisplacementDim);
    }

    return cache;
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int DisplacementDim>
std::vector<double> THMLocalAssembler<
    ShapeFunctionDisplacement, ShapeFunctionPressure, IntegrationMethod,
    DisplacementDim>::getDarcyVelocity() const
{
    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    std::vector<double> ip_q_values;
    auto cache_mat = MathLib::createZeroedMatrix<Eigen::Matrix<
        double, DisplacementDim, Eigen::Dynamic, Eigen::RowMajor>>(
        ip_q_values, DisplacementDim, n_integration_points);

    for (unsigned ip = 0; ip < n_integration_points; ++ip)
    {
        cache_mat.row(ip) = _ip_data[ip].q.transpose();
    }

    return ip_q_values;
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int DisplacementDim>
std::size_t THMLocalAssembler<
    ShapeFunctionDisplacement, ShapeFunctionPressure, IntegrationMethod,
    DisplacementDim>::setSigma(double const* values)
{
    auto const kelvin_vector_size =
        MathLib::KelvinVector::KelvinVectorDimensions<DisplacementDim>::value;
    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    std::vector<double> ip_sigma_values;
    auto sigma_values =
        Eigen::Map<Eigen::Matrix<double, kelvin_vector_size, Eigen::Dynamic,
                                 Eigen::ColMajor> const>(
            values, kelvin_vector_size, n_integration_points);

    for (unsigned ip = 0; ip < n_integration_points; ++ip)
    {
        _ip_data[ip].sigma =
            MathLib::KelvinVector::symmetricTensorToKelvinVector(
                sigma_values.col(ip));
    }

    return n_integration_points;
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int DisplacementDim>
std::vector<double> THMLocalAssembler<
    ShapeFunctionDisplacement, ShapeFunctionPressure, IntegrationMethod,
    DisplacementDim>::getSigma() const
{
    auto const kelvin_vector_size =
        MathLib::KelvinVector::KelvinVectorDimensions<DisplacementDim>::value;
    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    std::vector<double> ip_sigma_values;
    auto cache_mat = MathLib::createZeroedMatrix<Eigen::Matrix<
        double, Eigen::Dynamic, kelvin_vector_size, Eigen::RowMajor>>(
        ip_sigma_values, n_integration_points, kelvin_vector_size);

    for (unsigned ip = 0; ip < n_integration_points; ++ip)
    {
        auto const& sigma = _ip_data[ip].sigma_eff;
        cache_mat.row(ip) =
            MathLib::KelvinVector::kelvinVectorToSymmetricTensor(sigma);
    }

    return ip_sigma_values;
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int DisplacementDim>
std::vector<double> const& THMLocalAssembler<
    ShapeFunctionDisplacement, ShapeFunctionPressure, IntegrationMethod,
    DisplacementDim>::
    getIntPtSigma(const double /*t*/,
                  GlobalVector const& /*current_solution*/,
                  NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
                  std::vector<double>& cache) const
{
    static const int kelvin_vector_size =
        MathLib::KelvinVector::KelvinVectorDimensions<DisplacementDim>::value;
    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    cache.clear();
    auto cache_mat = MathLib::createZeroedMatrix<Eigen::Matrix<
        double, kelvin_vector_size, Eigen::Dynamic, Eigen::RowMajor>>(
        cache, kelvin_vector_size, n_integration_points);

    for (unsigned ip = 0; ip < n_integration_points; ++ip)
    {
        auto const& sigma = _ip_data[ip].sigma_eff;
        cache_mat.col(ip) =
            MathLib::KelvinVector::kelvinVectorToSymmetricTensor(sigma);
    }

    return cache;
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int DisplacementDim>
std::size_t THMLocalAssembler<
    ShapeFunctionDisplacement, ShapeFunctionPressure, IntegrationMethod,
    DisplacementDim>::setEpsilon(double const* values)
{
    auto const kelvin_vector_size =
        MathLib::KelvinVector::KelvinVectorDimensions<DisplacementDim>::value;
    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    auto epsilon_values =
        Eigen::Map<Eigen::Matrix<double, kelvin_vector_size, Eigen::Dynamic,
                                 Eigen::ColMajor> const>(
            values, kelvin_vector_size, n_integration_points);

    for (unsigned ip = 0; ip < n_integration_points; ++ip)
    {
        _ip_data[ip].eps = MathLib::KelvinVector::symmetricTensorToKelvinVector(
            epsilon_values.col(ip));
    }

    return n_integration_points;
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int DisplacementDim>
std::vector<double> THMLocalAssembler<
    ShapeFunctionDisplacement, ShapeFunctionPressure, IntegrationMethod,
    DisplacementDim>::getEpsilon() const
{
    auto const kelvin_vector_size =
        MathLib::KelvinVector::KelvinVectorDimensions<DisplacementDim>::value;
    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    std::vector<double> ip_epsilon_values;
    auto cache_mat = MathLib::createZeroedMatrix<Eigen::Matrix<
        double, Eigen::Dynamic, kelvin_vector_size, Eigen::RowMajor>>(
        ip_epsilon_values, n_integration_points, kelvin_vector_size);

    for (unsigned ip = 0; ip < n_integration_points; ++ip)
    {
        auto const& eps = _ip_data[ip].eps;
        cache_mat.row(ip) =
            MathLib::KelvinVector::kelvinVectorToSymmetricTensor(eps);
    }

    return ip_epsilon_values;
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int DisplacementDim>
std::vector<double> const& THMLocalAssembler<
    ShapeFunctionDisplacement, ShapeFunctionPressure, IntegrationMethod,
    DisplacementDim>::
    getIntPtEpsilon(const double /*t*/,
                    GlobalVector const& /*current_solution*/,
                    NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
                    std::vector<double>& cache) const
{
    auto const kelvin_vector_size =
        MathLib::KelvinVector::KelvinVectorDimensions<DisplacementDim>::value;
    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    cache.clear();
    auto cache_mat = MathLib::createZeroedMatrix<Eigen::Matrix<
        double, kelvin_vector_size, Eigen::Dynamic, Eigen::RowMajor>>(
        cache, kelvin_vector_size, n_integration_points);

    for (unsigned ip = 0; ip < n_integration_points; ++ip)
    {
        auto const& eps = _ip_data[ip].eps;
        cache_mat.col(ip) =
            MathLib::KelvinVector::kelvinVectorToSymmetricTensor(eps);
    }

    return cache;
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int DisplacementDim>
std::size_t THMLocalAssembler<
    ShapeFunctionDisplacement, ShapeFunctionPressure, IntegrationMethod,
    DisplacementDim>::setEpsilonMechanical(double const* values)
{
    auto const kelvin_vector_size =
        MathLib::KelvinVector::KelvinVectorDimensions<DisplacementDim>::value;
    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    auto epsilon_m_values =
        Eigen::Map<Eigen::Matrix<double, kelvin_vector_size, Eigen::Dynamic,
                                 Eigen::ColMajor> const>(
            values, kelvin_vector_size, n_integration_points);

    for (unsigned ip = 0; ip < n_integration_points; ++ip)
    {
        _ip_data[ip].eps_m =
            MathLib::KelvinVector::symmetricTensorToKelvinVector(
                epsilon_m_values.col(ip));
    }

    return n_integration_points;
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int DisplacementDim>
std::vector<double> THMLocalAssembler<
    ShapeFunctionDisplacement, ShapeFunctionPressure, IntegrationMethod,
    DisplacementDim>::getEpsilonMechanical() const
{
    auto const kelvin_vector_size =
        MathLib::KelvinVector::KelvinVectorDimensions<DisplacementDim>::value;
    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    std::vector<double> ip_epsilon_m_values;
    auto cache_mat = MathLib::createZeroedMatrix<Eigen::Matrix<
        double, Eigen::Dynamic, kelvin_vector_size, Eigen::RowMajor>>(
        ip_epsilon_m_values, n_integration_points, kelvin_vector_size);

    for (unsigned ip = 0; ip < n_integration_points; ++ip)
    {
        auto const& eps_m = _ip_data[ip].eps_m;
        cache_mat.row(ip) =
            MathLib::KelvinVector::kelvinVectorToSymmetricTensor(eps_m);
    }

    return ip_epsilon_m_values;
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int DisplacementDim>
void THMLocalAssembler<ShapeFunctionDisplacement,
                                        ShapeFunctionPressure,
                                        IntegrationMethod, DisplacementDim>::
    computeSecondaryVariableConcrete(double const /*t*/,
                                     std::vector<double> const& local_x)
{
    //auto p = Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
    //    pressure_size> const>(local_x.data() + pressure_index, pressure_size);

    //NumLib::interpolateToHigherOrderNodes<
    //    ShapeFunctionPressure, typename ShapeFunctionDisplacement::MeshElement,
    //    DisplacementDim>(_element, _is_axially_symmetric, p,
    //                     *_process_data.pressure_interpolated);

    //auto T = Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
    //    pressure_size> const>(local_x.data() + temperature_index, temperature_size);

    //NumLib::interpolateToHigherOrderNodes<
    //    ShapeFunctionPressure, typename ShapeFunctionDisplacement::MeshElement,
    //    DisplacementDim>(_element, _is_axially_symmetric, T,
    //                     *_process_data.temperature_interpolated);
}

}  // namespace THM
}  // namespace ProcessLib
