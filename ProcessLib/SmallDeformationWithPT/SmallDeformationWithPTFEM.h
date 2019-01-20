/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <memory>
#include <vector>

#include "MaterialLib/SolidModels/SelectSolidConstitutiveRelation.h"
#include "MathLib/KelvinVector.h"
#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "NumLib/Extrapolation/ExtrapolatableElement.h"
#include "NumLib/Fem/FiniteElement/TemplateIsoparametric.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "ProcessLib/Deformation/BMatrixPolicy.h"
#include "ProcessLib/Deformation/LinearBMatrix.h"
#include "ProcessLib/Parameter/Parameter.h"
#include "ProcessLib/Utils/InitShapeMatrices.h"

#include "LocalAssemblerInterface.h"
#include "SmallDeformationWithPTProcessData.h"

namespace ProcessLib
{
namespace SmallDeformationWithPT
{
template <typename BMatricesType, typename ShapeMatricesType,
          int DisplacementDim>
struct IntegrationPointData final
{
    explicit IntegrationPointData(
        MaterialLib::Solids::MechanicsBase<DisplacementDim> const&
            solid_material)
        : solid_material(solid_material),
          material_state_variables(
              solid_material.createMaterialStateVariables())
    {
    }

    // total stress
    typename BMatricesType::KelvinVectorType sigma, sigma_prev;

    // effective stress
    typename BMatricesType::KelvinVectorType eff_sigma, eff_sigma_prev;

    // total strain
    typename BMatricesType::KelvinVectorType eps, eps_prev;

    // mechanical strain
    typename BMatricesType::KelvinVectorType eps_m, eps_m_prev;

    MaterialLib::Solids::MechanicsBase<DisplacementDim> const& solid_material;
    std::unique_ptr<typename MaterialLib::Solids::MechanicsBase<
        DisplacementDim>::MaterialStateVariables>
        material_state_variables;
    // double solid_density;
    // double solid_density_prev;

    double integration_weight;
    typename ShapeMatricesType::NodalRowVectorType N;
    typename ShapeMatricesType::GlobalDimNodalMatrixType dNdx;

    void pushBackState()
    {
        eps_prev = eps;
        eps_m_prev = eps_m;
        sigma_prev = sigma;
        eff_sigma_prev = eff_sigma;
        // solid_density_prev = solid_density;
        material_state_variables->pushBackState();
    }

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

/// Used for the extrapolation of the integration point values. It is ordered
/// (and stored) by integration points.
template <typename ShapeMatrixType>
struct SecondaryData
{
    std::vector<ShapeMatrixType, Eigen::aligned_allocator<ShapeMatrixType>> N;
};

template <typename ShapeFunction, typename IntegrationMethod,
          int DisplacementDim>
class SmallDeformationWithPTLocalAssembler
    : public SmallDeformationWithPTLocalAssemblerInterface
{
public:
    using ShapeMatricesType =
        ShapeMatrixPolicyType<ShapeFunction, DisplacementDim>;

    // Types for displacement.
    // (Higher order elements = ShapeFunction).
    using ShapeMatrices = typename ShapeMatricesType::ShapeMatrices;
    using BMatricesType = BMatrixPolicyType<ShapeFunction, DisplacementDim>;

    using NodalForceVectorType = typename BMatricesType::NodalForceVectorType;
    using RhsVector = typename ShapeMatricesType::template VectorType<ShapeFunction::NPOINTS * DisplacementDim>;
    using JacobianMatrix = typename ShapeMatricesType::template MatrixType<
        ShapeFunction::NPOINTS * DisplacementDim,
        ShapeFunction::NPOINTS * DisplacementDim>;

    SmallDeformationWithPTLocalAssembler(SmallDeformationWithPTLocalAssembler const&) =
        delete;
    SmallDeformationWithPTLocalAssembler(SmallDeformationWithPTLocalAssembler&&) = delete;

    SmallDeformationWithPTLocalAssembler(
        MeshLib::Element const& e,
        std::size_t const /*local_matrix_size*/,
        bool const is_axially_symmetric,
        unsigned const integration_order,
        SmallDeformationWithPTProcessData<DisplacementDim>& process_data)
        : _process_data(process_data),
          _integration_method(integration_order),
          _element(e),
          _is_axially_symmetric(is_axially_symmetric)
    {
        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        _ip_data.reserve(n_integration_points);
        _secondary_data.N.resize(n_integration_points);

        auto const shape_matrices =
            initShapeMatrices<ShapeFunction, ShapeMatricesType,
                              IntegrationMethod, DisplacementDim>(
                e, is_axially_symmetric, _integration_method);

        auto& solid_material =
            MaterialLib::Solids::selectSolidConstitutiveRelation(
                _process_data.solid_materials,
                _process_data.material_ids,
                e.getID());

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            _ip_data.emplace_back(solid_material);
            auto& ip_data = _ip_data[ip];
            ip_data.integration_weight =
                _integration_method.getWeightedPoint(ip).getWeight() *
                shape_matrices[ip].integralMeasure * shape_matrices[ip].detJ;

            static const int kelvin_vector_size =
                MathLib::KelvinVector::KelvinVectorDimensions<
                    DisplacementDim>::value;
            ip_data.sigma.setZero(kelvin_vector_size);
            ip_data.sigma_prev.setZero(kelvin_vector_size);
            ip_data.eff_sigma.setZero(kelvin_vector_size);
            ip_data.eff_sigma_prev.setZero(kelvin_vector_size);
            ip_data.eps.setZero(kelvin_vector_size);
            ip_data.eps_prev.setZero(kelvin_vector_size);
            ip_data.eps_m.setZero(kelvin_vector_size);
            ip_data.eps_m_prev.setZero(kelvin_vector_size);

            // SpatialPosition x_position;
            // x_position.setElementID(_element.getID());
            // ip_data.solid_density =
            //     _process_data.reference_solid_density(0, x_position)[0];
            // ip_data.solid_density_prev = ip_data.solid_density;

            ip_data.N = shape_matrices[ip].N;
            ip_data.dNdx = shape_matrices[ip].dNdx;

            _secondary_data.N[ip] = shape_matrices[ip].N;
        }
    }

    /// Returns number of read integration points.
    std::size_t setIPDataInitialConditions(std::string const& name,
                                           double const* values,
                                           int const integration_order) override
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

    void assemble(double const /*t*/, std::vector<double> const& /*local_x*/,
                  std::vector<double>& /*local_M_data*/,
                  std::vector<double>& /*local_K_data*/,
                  std::vector<double>& /*local_rhs_data*/) override
    {
        OGS_FATAL(
            "SmallDeformationWithPTLocalAssembler: assembly without jacobian is not "
            "implemented.");
    }

    void assembleWithJacobian(double const t,
                              std::vector<double> const& local_x,
                              std::vector<double> const& /*local_xdot*/,
                              const double /*dxdot_dx*/, const double /*dx_dx*/,
                              std::vector<double>& /*local_M_data*/,
                              std::vector<double>& /*local_K_data*/,
                              std::vector<double>& local_rhs_data,
                              std::vector<double>& local_Jac_data) override
    {
        auto const local_matrix_size = local_x.size();
        assert(local_matrix_size == displacement_size);

        auto u = Eigen::Map<typename ShapeMatricesType::template VectorType<
            displacement_size> const>(local_x.data(), displacement_size);

        auto local_Jac = MathLib::createZeroedMatrix<JacobianMatrix>(
            local_Jac_data, local_matrix_size, local_matrix_size);

        auto local_rhs = MathLib::createZeroedVector<RhsVector>(
            local_rhs_data, local_matrix_size);

        // // set externaly given variables
        // Eigen::VectorXd nodal_T0(_element.getNumberOfNodes());
        // Eigen::VectorXd nodal_T1(_element.getNumberOfNodes());
        // Eigen::VectorXd nodal_p0(_element.getNumberOfNodes());
        // Eigen::VectorXd nodal_p1(_element.getNumberOfNodes());
        // for (unsigned i=0; i<_element.getNumberOfNodes(); i++)
        // {
        //     SpatialPosition x_position;
        //     x_position.setElementID(_element.getID());
        //     x_position.setNodeID(_element.getNode(i)->getID());
        //     nodal_T0[i] = _process_data.T0(t, x_position)[0];
        //     nodal_T1[i] = _process_data.T1(t, x_position)[0];
        //     nodal_p0[i] = _process_data.p0(t, x_position)[0];
        //     nodal_p1[i] = _process_data.p1(t, x_position)[0];
        // }

        double const& dt = _process_data.dt;

        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        SpatialPosition x_position;
        x_position.setElementID(_element.getID());

        using Invariants = MathLib::KelvinVector::Invariants<
            MathLib::KelvinVector::KelvinVectorDimensions<
                DisplacementDim>::value>;


        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            x_position.setIntegrationPoint(ip);
            auto const& w = _ip_data[ip].integration_weight;
            auto const& N = _ip_data[ip].N;
            auto const& dNdx = _ip_data[ip].dNdx;

            auto const x_coord =
                interpolateXCoordinate<ShapeFunction, ShapeMatricesType>(
                    _element, N);
            auto const& B = LinearBMatrix::computeBMatrix<
                DisplacementDim, ShapeFunction::NPOINTS,
                typename BMatricesType::BMatrixType>(dNdx, N, x_coord,
                                                     _is_axially_symmetric);

            auto& sigma = _ip_data[ip].sigma;
            auto const& sigma_prev = _ip_data[ip].sigma_prev;

            auto& eff_sigma = _ip_data[ip].eff_sigma;
            auto& eff_sigma_prev = _ip_data[ip].eff_sigma_prev;

            auto& eps = _ip_data[ip].eps;
            auto const& eps_prev = _ip_data[ip].eps_prev;

            auto& eps_m = _ip_data[ip].eps_m;
            auto const& eps_m_prev = _ip_data[ip].eps_m_prev;

            auto& state = _ip_data[ip].material_state_variables;

            auto const p0_ip = _process_data.p0(t, x_position)[0];
            auto const p1_ip = _process_data.p1(t, x_position)[0];
            auto const T0_ip = _process_data.T0(t, x_position)[0];
            auto const T1_ip = _process_data.T1(t, x_position)[0];
            // double const p0_ip = N.dot(nodal_p0);
            // double const p1_ip = N.dot(nodal_p1);
            // double const T0_ip = N.dot(nodal_T0);
            // double const T1_ip = N.dot(nodal_T1);

            auto const alpha =
                _process_data.solid_linear_thermal_expansion_coefficient(
                    t, x_position)[0];
            auto const biot = _process_data.biot_coefficient(t, x_position)[0];
            auto const rho_s = _process_data.solid_density(t, x_position)[0];
            auto const rho_f = _process_data.fluid_density(t, x_position)[0];
            auto const porosity = _process_data.porosity(t, x_position)[0];
            auto const rho_bulk = rho_s * (1. - porosity) + porosity * rho_f;
            auto const& b = _process_data.specific_body_force;


            //------------------------------------------------------
            // strain calculation
            //------------------------------------------------------
            eps.noalias() = B * u;

            // calculate thermally induced strain
            // assume isotropic thermal expansion
            double const linear_thermal_strain_increment = alpha * (T1_ip - T0_ip);

            eps_m.noalias() =
                eps_m_prev + eps - eps_prev -
                linear_thermal_strain_increment * Invariants::identity2;

            //------------------------------------------------------
            // stress, C calculation
            //------------------------------------------------------
            eff_sigma_prev.noalias() = sigma_prev + biot * p0_ip * Invariants::identity2;
            auto&& solution = _ip_data[ip].solid_material.integrateStress(
                t, x_position, dt, eps_m_prev, eps_m, eff_sigma_prev, *state, T1_ip);

            if (!solution)
                OGS_FATAL("Computation of local constitutive relation failed.");

            MathLib::KelvinVector::KelvinMatrixType<DisplacementDim> C;
            std::tie(eff_sigma, state, C) = std::move(*solution);
            sigma.noalias() = eff_sigma - biot * Invariants::identity2 * p1_ip;

            //------------------------------------------------------
            // residual, jacobian calculation
            //------------------------------------------------------
            local_Jac.noalias() += B.transpose() * C * B * w;

            typename ShapeMatricesType::template MatrixType<DisplacementDim,
                                                            displacement_size>
                N_u = ShapeMatricesType::template MatrixType<
                    DisplacementDim,
                    displacement_size>::Zero(DisplacementDim,
                                             displacement_size);

            for (int i = 0; i < DisplacementDim; ++i)
                N_u.template block<1, displacement_size / DisplacementDim>(
                       i, i * displacement_size / DisplacementDim)
                    .noalias() = N;

            //// calculate real density
            //// rho_s_{n+1} * (V_{n} + dV) = rho_s_n * V_n
            //// dV = 3 * alpha * dT * V_0
            //// rho_s_{n+1} = rho_s_n / (1 + 3 * alpha * dT )
            //// see reference solid density description for details.
            //auto& rho_s = _ip_data[ip].solid_density;
            //rho_s = _ip_data[ip].solid_density_prev /
            //        (1 + 3 * linear_thermal_strain_increment);

            local_rhs.noalias() -=
                (B.transpose() * sigma - N_u.transpose() * rho_bulk * b) * w;
        }
    }

    void preTimestepConcrete(std::vector<double> const& /*local_x*/,
                             double const /*t*/,
                             double const /*delta_t*/) override
    {
        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            _ip_data[ip].pushBackState();
        }
    }

    Eigen::Map<const Eigen::RowVectorXd> getShapeMatrix(
        const unsigned integration_point) const override
    {
        auto const& N = _secondary_data.N[integration_point];

        // assumes N is stored contiguously in memory
        return Eigen::Map<const Eigen::RowVectorXd>(N.data(), N.size());
    }

private:
    std::size_t setSigma(double const* values)
    {
        auto const kelvin_vector_size =
            MathLib::KelvinVector::KelvinVectorDimensions<
                DisplacementDim>::value;
        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

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

    // TODO (naumov) This method is same as getIntPtSigma but for arguments and
    // the ordering of the cache_mat.
    // There should be only one.
    std::vector<double> getSigma() const override
    {
        auto const kelvin_vector_size =
            MathLib::KelvinVector::KelvinVectorDimensions<
                DisplacementDim>::value;
        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        std::vector<double> ip_sigma_values;
        auto cache_mat = MathLib::createZeroedMatrix<Eigen::Matrix<
            double, Eigen::Dynamic, kelvin_vector_size, Eigen::RowMajor>>(
            ip_sigma_values, n_integration_points, kelvin_vector_size);

        for (unsigned ip = 0; ip < n_integration_points; ++ip)
        {
            auto const& sigma = _ip_data[ip].sigma;
            cache_mat.row(ip) =
                MathLib::KelvinVector::kelvinVectorToSymmetricTensor(sigma);
        }

        return ip_sigma_values;
    }

    std::vector<double> const& getIntPtSigma(
        const double /*t*/,
        GlobalVector const& /*current_solution*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        static const int kelvin_vector_size =
            MathLib::KelvinVector::KelvinVectorDimensions<
                DisplacementDim>::value;
        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        cache.clear();
        auto cache_mat = MathLib::createZeroedMatrix<Eigen::Matrix<
            double, kelvin_vector_size, Eigen::Dynamic, Eigen::RowMajor>>(
            cache, kelvin_vector_size, n_integration_points);

        for (unsigned ip = 0; ip < n_integration_points; ++ip)
        {
            auto const& sigma = _ip_data[ip].sigma;
            cache_mat.col(ip) =
                MathLib::KelvinVector::kelvinVectorToSymmetricTensor(sigma);
        }

        return cache;
    }

    std::size_t setEpsilon(double const* values)
    {
        auto const kelvin_vector_size =
            MathLib::KelvinVector::KelvinVectorDimensions<
                DisplacementDim>::value;
        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        auto epsilon_values =
            Eigen::Map<Eigen::Matrix<double, kelvin_vector_size, Eigen::Dynamic,
                                     Eigen::ColMajor> const>(
                values, kelvin_vector_size, n_integration_points);

        for (unsigned ip = 0; ip < n_integration_points; ++ip)
        {
            _ip_data[ip].eps =
                MathLib::KelvinVector::symmetricTensorToKelvinVector(
                    epsilon_values.col(ip));
        }

        return n_integration_points;
    }

    std::vector<double> getEpsilon() const override
    {
        auto const kelvin_vector_size =
            MathLib::KelvinVector::KelvinVectorDimensions<
                DisplacementDim>::value;
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

    virtual std::vector<double> const& getIntPtEpsilon(
        const double /*t*/,
        GlobalVector const& /*current_solution*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        auto const kelvin_vector_size =
            MathLib::KelvinVector::KelvinVectorDimensions<
                DisplacementDim>::value;
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

    std::size_t setEpsilonMechanical(double const* values)
    {
        auto const kelvin_vector_size =
            MathLib::KelvinVector::KelvinVectorDimensions<
                DisplacementDim>::value;
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

    std::vector<double> getEpsilonMechanical() const override
    {
        auto const kelvin_vector_size =
            MathLib::KelvinVector::KelvinVectorDimensions<
                DisplacementDim>::value;
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

    SmallDeformationWithPTProcessData<DisplacementDim>& _process_data;

    std::vector<
        IntegrationPointData<BMatricesType, ShapeMatricesType, DisplacementDim>,
        Eigen::aligned_allocator<IntegrationPointData<
            BMatricesType, ShapeMatricesType, DisplacementDim>>>
        _ip_data;

    IntegrationMethod _integration_method;
    MeshLib::Element const& _element;
    SecondaryData<typename ShapeMatrices::ShapeType> _secondary_data;
    bool const _is_axially_symmetric;

    static const int displacement_index = ShapeFunction::NPOINTS;
    static const int displacement_size =
        ShapeFunction::NPOINTS * DisplacementDim;
};

}  // namespace SmallDeformationWithPT
}  // namespace ProcessLib
