
#pragma once

#include <vector>

#include "NumLib/Fem/ShapeMatrixPolicy.h"

#include "ProcessLib/LIE/Common/HMatrixUtils.h"
#include "ProcessLib/LIE/TH/THProcessData.h"

#include "IntegrationPointDataFracture.h"
#include "THLocalAssemblerInterface.h"

namespace ProcessLib
{
namespace LIE
{
namespace TH
{
template <typename ShapeFunctionPressure,
          typename IntegrationMethod,
          int GlobalDim>
class THLocalAssemblerFracture : public THLocalAssemblerInterface
{
public:
    THLocalAssemblerFracture(THLocalAssemblerFracture const&) = delete;
    THLocalAssemblerFracture(THLocalAssemblerFracture&&) = delete;

    THLocalAssemblerFracture(
        MeshLib::Element const& e,
        std::size_t const n_variables,
        std::size_t const local_matrix_size,
        std::vector<unsigned> const& dofIndex_to_localIndex,
        bool const is_axially_symmetric,
        unsigned const integration_order,
        THProcessData<GlobalDim>& process_data);

    void preTimestepConcrete(std::vector<double> const& /*local_x*/,
                             double const /*t*/,
                             double const /*delta_t*/) override
    {
        for (auto& data : _ip_data)
        {
            data.pushBackState();
        }
    }

    void computeSecondaryVariableConcreteWithVector(
        const double t, Eigen::VectorXd const& local_x) override;

    Eigen::Map<const Eigen::RowVectorXd> getShapeMatrix(
        const unsigned integration_point) const override
    {
        auto const& N = _ip_data[integration_point].N_p;

        // assumes N is stored contiguously in memory
        return Eigen::Map<const Eigen::RowVectorXd>(N.data(), N.size());
    }

private:
    void assembleWithJacobianConcrete(double const t,
                                      Eigen::VectorXd const& local_x,
                                      Eigen::VectorXd const& local_xdot,
                                      Eigen::VectorXd& local_b,
                                      Eigen::MatrixXd& local_J) override;

    void assembleBlockMatricesWithJacobian(
        double const t,
        Eigen::Ref<const Eigen::VectorXd> const& p,
        Eigen::Ref<const Eigen::VectorXd> const& p_dot,
        Eigen::Ref<const Eigen::VectorXd> const& T,
        Eigen::Ref<const Eigen::VectorXd> const& T_dot,
        Eigen::Ref<Eigen::VectorXd>
            rhs_p,
        Eigen::Ref<Eigen::VectorXd>
            rhs_T,
        Eigen::Ref<Eigen::MatrixXd>
            J_pp,
        Eigen::Ref<Eigen::MatrixXd>
            J_pT,
        Eigen::Ref<Eigen::MatrixXd>
            J_TT,
        Eigen::Ref<Eigen::MatrixXd>
            J_Tp);

    // Types for pressure.
    using ShapeMatricesTypePressure =
        ShapeMatrixPolicyType<ShapeFunctionPressure, GlobalDim>;

    // Types for the integration point data
    using IntegrationPointDataType =
        IntegrationPointDataFracture<ShapeMatricesTypePressure,
                                     GlobalDim>;

    THProcessData<GlobalDim>& _process_data;
    std::vector<FractureProperty*> _fracture_props;
    std::vector<JunctionProperty*> _junction_props;
    std::unordered_map<int, int> _fracID_to_local;
    FractureProperty const* _fracture_property = nullptr;

    std::vector<IntegrationPointDataType,
                Eigen::aligned_allocator<IntegrationPointDataType>>
        _ip_data;

    Eigen::Vector3d _e_center_coords;

    static constexpr int pressure_index = 0;
    static constexpr int pressure_size = ShapeFunctionPressure::NPOINTS;
    static constexpr int temperature_index = pressure_index + pressure_size;
    static constexpr int temperature_size = ShapeFunctionPressure::NPOINTS;
};

}  // namespace TH
}  // namespace LIE
}  // namespace ProcessLib

#include "THLocalAssemblerFracture-impl.h"
