
#pragma once

#include <memory>
#include <vector>

#include "MathLib/KelvinVector.h"
#include "MaterialLib/SolidModels/MechanicsBase.h"

namespace ProcessLib
{
namespace LIE
{
namespace TH
{
template <typename ShapeMatrixTypePressure, unsigned GlobalDim,
          unsigned NPoints>
struct IntegrationPointDataMatrix final
{
    explicit IntegrationPointDataMatrix()
    {
        q.setZero();
    }

    typename ShapeMatrixTypePressure::NodalRowVectorType N_p;
    typename ShapeMatrixTypePressure::GlobalDimNodalMatrixType dNdx_p;

    double integration_weight;

    Eigen::Matrix<double, GlobalDim, 1> q;
    Eigen::Vector3d darcy_velocity;

    void pushBackState()
    {
    }

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

}  // namespace TH
}  // namespace LIE
}  // namespace ProcessLib
