
#pragma once

#include <Eigen/Eigen>

#include "MaterialLib/FractureModels/FractureModelBase.h"

namespace ProcessLib
{
namespace LIE
{
namespace TH
{
template <typename ShapeMatrixTypePressure, unsigned GlobalDim>
struct IntegrationPointDataFracture final
{
    explicit IntegrationPointDataFracture()
    {
    }

    typename ShapeMatrixTypePressure::NodalRowVectorType N_p;
    typename ShapeMatrixTypePressure::GlobalDimNodalMatrixType dNdx_p;

    double aperture = 0.0;
    double aperture0 = 0.0;
    double permeability = 0.0;
    Eigen::Matrix<double, GlobalDim, 1> q;

    std::unique_ptr<
        typename MaterialLib::Fracture::Permeability::PermeabilityState>
        permeability_state;

    double integration_weight;

    Eigen::Vector3d darcy_velocity;

    void pushBackState()
    {
    }

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

}  // namespace TH
}  // namespace LIE
}  // namespace ProcessLib
