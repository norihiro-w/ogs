
#pragma once

#include <memory>
#include <vector>

#include "MathLib/KelvinVector.h"
#include "MaterialLib/SolidModels/MechanicsBase.h"

namespace ProcessLib
{
namespace LIE
{
namespace ThermoHydroMechanics
{
template <typename BMatricesType, typename ShapeMatrixTypeDisplacement,
          typename ShapeMatrixTypePressure, unsigned GlobalDim,
          unsigned NPoints>
struct IntegrationPointDataMatrix final
{
    explicit IntegrationPointDataMatrix(
        MaterialLib::Solids::MechanicsBase<GlobalDim>& solid_material)
        : solid_material(solid_material),
          material_state_variables(
              solid_material.createMaterialStateVariables())
    {
        static const int kelvin_vector_size =
            MathLib::KelvinVector::KelvinVectorDimensions<
                GlobalDim>::value;

        sigma.setZero(kelvin_vector_size);
        sigma_prev.resize(kelvin_vector_size);
        sigma_eff.setZero(kelvin_vector_size);
        sigma_eff_prev.resize(kelvin_vector_size);
        eps.setZero(kelvin_vector_size);
        eps_prev.resize(kelvin_vector_size);
        eps_m.setZero(kelvin_vector_size);
        eps_m_prev.resize(kelvin_vector_size);
        q.setZero();
    }

    typename ShapeMatrixTypeDisplacement::NodalRowVectorType N_u;
    typename ShapeMatrixTypeDisplacement::GlobalDimNodalMatrixType dNdx_u;
    typename ShapeMatrixTypeDisplacement::template MatrixType<
        GlobalDim, NPoints * GlobalDim>
        H_u;
    typename BMatricesType::KelvinVectorType sigma, sigma_prev;
    typename BMatricesType::KelvinVectorType sigma_eff, sigma_eff_prev;
    typename BMatricesType::KelvinVectorType eps, eps_prev;
    typename BMatricesType::KelvinVectorType eps_m, eps_m_prev;
    Eigen::Matrix<double, GlobalDim, 1> q;

    typename ShapeMatrixTypePressure::NodalRowVectorType N_p;
    typename ShapeMatrixTypePressure::GlobalDimNodalMatrixType dNdx_p;

    MaterialLib::Solids::MechanicsBase<GlobalDim>& solid_material;
    std::unique_ptr<typename MaterialLib::Solids::MechanicsBase<
        GlobalDim>::MaterialStateVariables>
        material_state_variables;

    typename BMatricesType::KelvinMatrixType C;
    double integration_weight;

    Eigen::Vector3d darcy_velocity;

    void pushBackState()
    {
        eps_prev = eps;
        eps_m_prev = eps_m;
        sigma_prev = sigma;
        sigma_eff_prev = sigma_eff;
        material_state_variables->pushBackState();
    }

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

}  // namespace ThermoHydroMechanics
}  // namespace LIE
}  // namespace ProcessLib
