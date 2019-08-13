/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "HydroMechanicsLocalAssemblerFracture.h"

#include "MaterialLib/FractureModels/FractureIdentity2.h"

#include "ProcessLib/Utils/InitShapeMatrices.h"

#include "ProcessLib/LIE/Common/LevelSetFunction.h"

namespace ProcessLib
{
namespace LIE
{
namespace HydroMechanics
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

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int GlobalDim>
HydroMechanicsLocalAssemblerFracture<ShapeFunctionDisplacement,
                                     ShapeFunctionPressure, IntegrationMethod,
                                     GlobalDim>::
    HydroMechanicsLocalAssemblerFracture(
        MeshLib::Element const& e,
        std::size_t const n_variables,
        std::size_t const /*local_matrix_size*/,
        std::vector<unsigned> const& dofIndex_to_localIndex,
        bool const is_axially_symmetric,
        unsigned const integration_order,
        HydroMechanicsProcessData<GlobalDim>& process_data)
    : HydroMechanicsLocalAssemblerInterface(
          e, is_axially_symmetric,
          n_variables * ShapeFunctionDisplacement::NPOINTS * GlobalDim +
              ShapeFunctionPressure::NPOINTS,
          dofIndex_to_localIndex),
      _process_data(process_data),
      _e_center_coords(e.getCenterOfGravity().getCoords())
{
    assert(e.getDimension() == GlobalDim - 1);

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

    // Get element nodes for aperture0 interpolation from nodes to integration
    // point. The aperture0 parameter is time-independent.
    typename ShapeMatricesTypeDisplacement::NodalVectorType
        aperture0_node_values =
            _fracture_property->aperture0.getNodalValuesOnElement(
            e, /*time independent*/ 0);

    ParameterLib::SpatialPosition x_position;
    x_position.setElementID(e.getID());
    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        x_position.setIntegrationPoint(ip);

        _ip_data.emplace_back(*_process_data.fracture_model);
        auto const& sm_u = shape_matrices_u[ip];
        auto const& sm_p = shape_matrices_p[ip];
        auto& ip_data = _ip_data[ip];
        ip_data.integration_weight =
            sm_u.detJ * sm_u.integralMeasure *
            integration_method.getWeightedPoint(ip).getWeight();

        ip_data.H_u.setZero(GlobalDim,
                            ShapeFunctionDisplacement::NPOINTS * GlobalDim);
        computeHMatrix<
            GlobalDim, ShapeFunctionDisplacement::NPOINTS,
            typename ShapeMatricesTypeDisplacement::NodalRowVectorType,
            HMatrixType>(sm_u.N, ip_data.H_u);
        ip_data.N_p = sm_p.N;
        ip_data.dNdx_p = sm_p.dNdx;

        // Initialize current time step values
        ip_data.w.setZero(GlobalDim);
        ip_data.sigma_eff.setZero(GlobalDim);

        // Previous time step values are not initialized and are set later.
        ip_data.w_prev.resize(GlobalDim);
        ip_data.sigma_eff_prev.resize(GlobalDim);

        ip_data.C.resize(GlobalDim, GlobalDim);

        ip_data.aperture0 = aperture0_node_values.dot(sm_u.N);
        ip_data.aperture = ip_data.aperture0;

        ip_data.permeability_state =
            _fracture_property->permeability_model->getNewState();

        auto const initial_effective_stress = (*_fracture_property->initial_fracture_effective_stress)(0, x_position);
        for (int i = 0; i < GlobalDim; i++)
        {
            ip_data.sigma_eff[i] = initial_effective_stress[i];
            ip_data.sigma_eff_prev[i] = initial_effective_stress[i];
        }
    }
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int GlobalDim>
void HydroMechanicsLocalAssemblerFracture<
    ShapeFunctionDisplacement, ShapeFunctionPressure, IntegrationMethod,
    GlobalDim>::assembleWithJacobianConcrete(double const t,
                                             Eigen::VectorXd const& local_x,
                                             Eigen::VectorXd const& local_xdot,
                                             Eigen::VectorXd& local_b,
                                             Eigen::MatrixXd& local_J)
{
    auto const pressure_index = pressure_index_;
    auto const pressure_size = pressure_size_;
    auto const displacement_jump_index = displacement_jump_index_;
    auto const displacement_jump_size = displacement_jump_size_;
    auto const n_fractures = _fracture_props.size();
    auto const n_junctions = _junction_props.size();
    auto const n_enrich_var = n_fractures + n_junctions;

    auto const local_p = local_x.segment<pressure_size>(pressure_index);
    auto const local_p_dot = local_xdot.segment<pressure_size>(pressure_index);
    std::vector<Eigen::VectorXd> vec_local_g;
    std::vector<Eigen::VectorXd> vec_local_g_dot;
    for (unsigned i = 0; i < n_enrich_var; i++)
    {
        auto sub = local_x.segment<displacement_jump_size>(displacement_jump_index + displacement_jump_size * i);
        vec_local_g.push_back(sub);
        auto sub_dot = local_xdot.segment<displacement_jump_size>(displacement_jump_index + displacement_jump_size * i);
        vec_local_g_dot.push_back(sub_dot);
    }

    //
    using BlockVectorTypeU =
        typename Eigen::VectorXd::FixedSegmentReturnType<displacement_jump_size>::Type;
    using BlockMatrixTypePU =
        Eigen::Block<Eigen::MatrixXd, pressure_size, displacement_jump_size>;
    using BlockMatrixTypeUP =
        Eigen::Block<Eigen::MatrixXd, displacement_jump_size, pressure_size>;
    using BlockMatrixTypeUU =
        Eigen::Block<Eigen::MatrixXd, displacement_jump_size, displacement_jump_size>;

    auto local_b_p = local_b.segment<pressure_size>(pressure_index);
    std::vector<BlockVectorTypeU> vec_local_b_g;
    for (unsigned i = 0; i < n_enrich_var; i++)
    {
        vec_local_b_g.push_back(
            local_b.segment<displacement_jump_size>(displacement_jump_index + displacement_jump_size*i));
    }
    auto local_J_pp = local_J.block<pressure_size,pressure_size>(pressure_index, pressure_index);
    std::vector<BlockMatrixTypePU> vec_local_J_pg;
    for (unsigned i = 0; i < n_enrich_var; i++)
    {
        auto sub_pg = local_J.block<pressure_size, displacement_jump_size>(
            pressure_index, displacement_jump_index + displacement_jump_size * i);
        vec_local_J_pg.push_back(sub_pg);
    }
    std::vector<BlockMatrixTypeUP> vec_local_J_gp;
    for (unsigned i = 0; i < n_enrich_var; i++)
    {
        auto sub_gp = local_J.block<displacement_jump_size, pressure_size>(
            displacement_jump_index + displacement_jump_size * i, pressure_index);
        vec_local_J_gp.push_back(sub_gp);
    }
    std::vector<std::vector<BlockMatrixTypeUU>> vec_local_J_gg(n_enrich_var);
    for (unsigned i = 0; i < n_enrich_var; i++)
    {
        for (unsigned j = 0; j < n_enrich_var; j++)
        {
            auto sub_gg = local_J.block<displacement_jump_size, displacement_jump_size>(
                displacement_jump_index + displacement_jump_size * i, displacement_jump_index + displacement_jump_size * j);
            vec_local_J_gg[i].push_back(sub_gg);
        }
    }

    // construct nodal w
    Eigen::VectorXd local_w(displacement_jump_size), local_w_dot(displacement_jump_size);
    local_w.setZero();
    local_w_dot.setZero();

    std::vector<double> const levelsets(duGlobalEnrichments(
        _fracture_property->fracture_id, _fracture_props, _junction_props,
        _fracID_to_local, _e_center_coords));
    for (unsigned i = 0; i < n_enrich_var; i++)
    {
        local_w.noalias() += levelsets[i] * vec_local_g[i];
        local_w_dot.noalias() += levelsets[i] * vec_local_g_dot[i];
    }

    Eigen::VectorXd local_b_w(displacement_jump_size);
    Eigen::MatrixXd local_J_pw(pressure_size, displacement_jump_size);
    Eigen::MatrixXd local_J_ww(displacement_jump_size, displacement_jump_size);
    Eigen::MatrixXd local_J_wp(displacement_jump_size, pressure_size);
    local_b_w.setZero();
    local_J_pw.setZero();
    local_J_ww.setZero();
    local_J_wp.setZero();

    assembleBlockMatricesWithJacobian(t, local_p, local_p_dot, local_w, local_w_dot,
                                      local_b_p, local_b_w, local_J_pp, local_J_pw,
                                      local_J_ww, local_J_wp);

    for (unsigned i = 0; i < n_enrich_var; i++)
    {
        vec_local_b_g[i].noalias() =  levelsets[i] * local_b_w;
        vec_local_J_pg[i].noalias() =  levelsets[i] * local_J_pw;
        vec_local_J_gp[i].noalias() =  levelsets[i] * local_J_wp;
        for (unsigned j = 0; j < n_enrich_var; j++)
            vec_local_J_gg[i][j].noalias() = levelsets[i] * levelsets[j] * local_J_ww;
    }

}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int GlobalDim>
void HydroMechanicsLocalAssemblerFracture<ShapeFunctionDisplacement,
                                          ShapeFunctionPressure,
                                          IntegrationMethod, GlobalDim>::
    assembleBlockMatricesWithJacobian(
        double const t, Eigen::Ref<const Eigen::VectorXd> const& p,
        Eigen::Ref<const Eigen::VectorXd> const& p_dot,
        Eigen::Ref<const Eigen::VectorXd> const& g,
        Eigen::Ref<const Eigen::VectorXd> const& g_dot,
        Eigen::Ref<Eigen::VectorXd> rhs_p, Eigen::Ref<Eigen::VectorXd> rhs_g,
        Eigen::Ref<Eigen::MatrixXd> J_pp, Eigen::Ref<Eigen::MatrixXd> J_pg,
        Eigen::Ref<Eigen::MatrixXd> J_gg, Eigen::Ref<Eigen::MatrixXd> J_gp)
{
    auto const pressure_size = pressure_size_;
    auto const displacement_jump_size = displacement_jump_size_;
    auto const& R = _fracture_property->R;
    double const& dt = _process_data.dt;

    // the index of a normal (normal to a fracture plane) component
    // in a displacement vector
    auto constexpr index_normal = GlobalDim - 1;

    typename ShapeMatricesTypePressure::NodalMatrixType laplace_p =
        ShapeMatricesTypePressure::NodalMatrixType::Zero(pressure_size,
                                                         pressure_size);

    typename ShapeMatricesTypePressure::NodalMatrixType storage_p =
        ShapeMatricesTypePressure::NodalMatrixType::Zero(pressure_size,
                                                         pressure_size);

    typename ShapeMatricesTypeDisplacement::template MatrixType<
        displacement_jump_size, pressure_size>
        Kgp = ShapeMatricesTypeDisplacement::template MatrixType<
            displacement_jump_size, pressure_size>::Zero(displacement_jump_size,
                                                    pressure_size);

    using GlobalDimMatrix = Eigen::Matrix<double, GlobalDim, GlobalDim>;
    using GlobalDimVector = Eigen::Matrix<double, GlobalDim, 1>;

    auto const& gravity_vec = _process_data.specific_body_force;

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
        auto const& H_g = ip_data.H_u;
        auto const& identity2 =
            MaterialLib::Fracture::FractureIdentity2<GlobalDim>::value;

        auto& mat = ip_data.fracture_material;
        auto& effective_stress = ip_data.sigma_eff;
        auto const& effective_stress_prev = ip_data.sigma_eff_prev;
        auto& w = ip_data.w;
        auto const& w_prev = ip_data.w_prev;
        auto& C = ip_data.C;
        auto& state = *ip_data.material_state_variables;
        auto& b_m = ip_data.aperture;

        double const S = (*_fracture_property->specific_storage)(t, x_position)[0];
        double const mu = _process_data.fluid_viscosity(t, x_position)[0];
        auto const alpha = (*_fracture_property->biot_coefficient)(t, x_position)[0];
        auto const rho_fr = _process_data.fluid_density(t, x_position)[0];

        // displacement jumps in local coordinates
        w.noalias() = R * H_g * g;

        // aperture
        b_m = ip_data.aperture0 + w[index_normal];
        if (b_m < 0.0)
        {
            DBUG(
                "Element %d, gp %d: Fracture aperture is %g, but it must be "
                "non-negative. Setting it to zero.",
                _element.getID(), ip, b_m);
            b_m = 0;
        }

        auto const initial_effective_stress =
            (*_fracture_property->initial_fracture_effective_stress)(0, x_position);

        Eigen::Map<typename HMatricesType::ForceVectorType const> const stress0(
            initial_effective_stress.data(), initial_effective_stress.size());

        // local C, local stress
        mat.computeConstitutiveRelation(
            t, x_position, ip_data.aperture0, stress0, w_prev, w,
            effective_stress_prev, effective_stress, C, state);

        auto& permeability = ip_data.permeability;
        permeability = _fracture_property->permeability_model->permeability(
            ip_data.permeability_state.get(), ip_data.aperture0, b_m);

        GlobalDimMatrix const k =
            createRotatedTensor<GlobalDim>(R, permeability);

        // derivative of permeability respect to aperture
        double const local_dk_db =
            _fracture_property->permeability_model->dpermeability_daperture(
                ip_data.permeability_state.get(), ip_data.aperture0, b_m);
        GlobalDimMatrix const dk_db =
            createRotatedTensor<GlobalDim>(R, local_dk_db);

        //
        // displacement equation, displacement jump part
        //
        rhs_g.noalias() -=
            H_g.transpose() * R.transpose() * effective_stress * ip_w;
        J_gg.noalias() += H_g.transpose() * R.transpose() * C * R * H_g * ip_w;

        //
        // displacement equation, pressure part
        //
        Kgp.noalias() +=
            H_g.transpose() * R.transpose() * alpha * identity2 * N_p * ip_w;

        //
        // pressure equation, pressure part.
        //
        storage_p.noalias() += N_p.transpose() * b_m * S * N_p * ip_w;
        laplace_p.noalias() +=
            dNdx_p.transpose() * b_m * k / mu * dNdx_p * ip_w;
        rhs_p.noalias() +=
            dNdx_p.transpose() * b_m * k / mu * rho_fr * gravity_vec * ip_w;

        //
        // pressure equation, displacement jump part.
        //
        GlobalDimVector const grad_head_over_mu =
            (dNdx_p * p + rho_fr * gravity_vec) / mu;
        Eigen::Matrix<double, 1, displacement_jump_size> const mT_R_Hg =
            identity2.transpose() * R * H_g;
        // velocity in global coordinates
        ip_data.darcy_velocity.head(GlobalDim).noalias() =
            -R.transpose() * k * grad_head_over_mu;
        J_pg.noalias() += N_p.transpose() * S * N_p * p_dot * mT_R_Hg * ip_w;
        J_pg.noalias() +=
            dNdx_p.transpose() * k * grad_head_over_mu * mT_R_Hg * ip_w;
        J_pg.noalias() += dNdx_p.transpose() * b_m * dk_db * grad_head_over_mu *
                          mT_R_Hg * ip_w;
    }

    // displacement equation, pressure part
    J_gp.noalias() -= Kgp;

    // pressure equation, pressure part.
    J_pp.noalias() += laplace_p + storage_p / dt;

    // pressure equation, displacement jump part.
    J_pg.noalias() += Kgp.transpose() / dt;

    // pressure equation
    rhs_p.noalias() -=
        laplace_p * p + storage_p * p_dot + Kgp.transpose() * g_dot;

    // displacement equation
    rhs_g.noalias() -= -Kgp * p;
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int GlobalDim>
void HydroMechanicsLocalAssemblerFracture<ShapeFunctionDisplacement,
                                          ShapeFunctionPressure,
                                          IntegrationMethod, GlobalDim>::
    computeSecondaryVariableConcreteWithVector(const double t,
                                               Eigen::VectorXd const& local_x)
{
    auto const displacement_jump_index = displacement_jump_index_;
    auto const displacement_jump_size = displacement_jump_size_;
    auto const n_fractures = _fracture_props.size();
    auto const n_junctions = _junction_props.size();
    auto const n_enrich_var = n_fractures + n_junctions;

    auto const& R = _fracture_property->R;

    // the index of a normal (normal to a fracture plane) component
    // in a displacement vector
    auto constexpr index_normal = GlobalDim - 1;

    ParameterLib::SpatialPosition x_position;
    x_position.setElementID(_element.getID());

    unsigned const n_integration_points = _ip_data.size();

    std::vector<Eigen::VectorXd> vec_nodal_g;
    for (unsigned i = 0; i < n_enrich_var; i++)
    {
        auto sub = local_x.segment<displacement_jump_size>(displacement_jump_index + displacement_jump_size * i);
        vec_nodal_g.push_back(sub);
    }

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        x_position.setIntegrationPoint(ip);

        auto& ip_data = _ip_data[ip];
        auto const& H_g = ip_data.H_u;
        auto& mat = ip_data.fracture_material;
        auto& effective_stress = ip_data.sigma_eff;
        auto const& effective_stress_prev = ip_data.sigma_eff_prev;
        auto& w = ip_data.w;
        auto const& w_prev = ip_data.w_prev;
        auto& C = ip_data.C;
        auto& state = *ip_data.material_state_variables;
        auto& b_m = ip_data.aperture;
        auto const& N = ip_data.N_p;

        Eigen::Vector3d const ip_physical_coords(
            computePhysicalCoordinates(_element, N).getCoords());
        std::vector<double> const levelsets(duGlobalEnrichments(
            _fracture_property->fracture_id, _fracture_props, _junction_props,
            _fracID_to_local, ip_physical_coords));

        Eigen::VectorXd nodal_gap(displacement_jump_size);
        nodal_gap.setZero();
        for (unsigned i = 0; i < n_enrich_var; i++)
            nodal_gap.noalias() += levelsets[i] * vec_nodal_g[i];

        // displacement jumps in local coordinates
        w.noalias() = R * H_g * nodal_gap;

        // aperture
        b_m = ip_data.aperture0 + w[index_normal];
        if (b_m < 0.0)
        {
            DBUG(
                "Element %d, gp %d: Fracture aperture is %g, but it is "
                "expected to be non-negative.",
                _element.getID(), ip, b_m);
        }

        auto const initial_effective_stress =
            (*_fracture_property->initial_fracture_effective_stress)(0, x_position);

        Eigen::Map<typename HMatricesType::ForceVectorType const> const stress0(
            initial_effective_stress.data(), initial_effective_stress.size());

        // local C, local stress
        mat.computeConstitutiveRelation(
            t, x_position, ip_data.aperture0, stress0, w_prev, w,
            effective_stress_prev, effective_stress, C, state);
    }

    double ele_b = 0;
    double ele_k = 0;
    typename HMatricesType::ForceVectorType ele_sigma_eff =
        HMatricesType::ForceVectorType::Zero(GlobalDim);
    typename HMatricesType::ForceVectorType ele_w =
        HMatricesType::ForceVectorType::Zero(GlobalDim);
    double ele_Fs = -std::numeric_limits<double>::max();
    Eigen::Vector3d ele_velocity = Eigen::Vector3d::Zero();
    for (auto const& ip : _ip_data)
    {
        ele_b += ip.aperture;
        ele_k += ip.permeability;
        ele_w += ip.w;
        ele_sigma_eff += ip.sigma_eff;
        ele_Fs = std::max(
            ele_Fs, ip.material_state_variables->getShearYieldFunctionValue());
        ele_velocity += ip.darcy_velocity;
    }
    ele_b /= static_cast<double>(n_integration_points);
    ele_k /= static_cast<double>(n_integration_points);
    ele_w /= static_cast<double>(n_integration_points);
    ele_sigma_eff /= static_cast<double>(n_integration_points);
    ele_velocity /= static_cast<double>(n_integration_points);
    auto const element_id = _element.getID();
    (*_process_data.mesh_prop_b)[element_id] = ele_b;
    (*_process_data.mesh_prop_k_f)[element_id] = ele_k;
    (*_process_data.mesh_prop_w_n)[element_id] = ele_w[index_normal];
    (*_process_data.mesh_prop_w_s)[element_id] = ele_w[0];
    (*_process_data.mesh_prop_fracture_stress_normal)[element_id] =
        ele_sigma_eff[index_normal];
    (*_process_data.mesh_prop_fracture_stress_shear)[element_id] =
        ele_sigma_eff[0];
    (*_process_data.mesh_prop_fracture_shear_failure)[element_id] = ele_Fs;

    if (GlobalDim == 3)
    {
        (*_process_data.mesh_prop_w_s2)[element_id] = ele_w[1];
        (*_process_data.mesh_prop_fracture_stress_shear2)[element_id] =
            ele_sigma_eff[1];
    }

    for (unsigned i = 0; i < 3; i++)
    {
        (*_process_data.mesh_prop_velocity)[element_id * 3 + i] =
            ele_velocity[i];
    }
}

}  // namespace HydroMechanics
}  // namespace LIE
}  // namespace ProcessLib
