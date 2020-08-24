
#pragma once

#include "ProcessLib/Process.h"

#include "ThermoHydroMechanicsProcessData.h"
#include "LocalAssembler/ThermoHydroMechanicsLocalAssemblerInterface.h"

namespace ProcessLib
{
namespace LIE
{
namespace ThermoHydroMechanics
{
class ThermoHydroMechanicsLocalAssemblerInterface;

template <int GlobalDim>
class ThermoHydroMechanicsProcess final : public Process
{
    static_assert(GlobalDim == 2 || GlobalDim == 3,
                  "Currently LIE::ThermoHydroMechanicsProcess "
                  "supports only 2D or 3D.");

public:
    ThermoHydroMechanicsProcess(
        std::string name,
        MeshLib::Mesh& mesh,
        std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&&
            jacobian_assembler,
        std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const&
            parameters,
        unsigned const integration_order,
        std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>&&
            process_variables,
        ThermoHydroMechanicsProcessData<GlobalDim>&& process_data,
        SecondaryVariableCollection&& secondary_variables,
        NumLib::NamedFunctionCaller&& named_function_caller,
        bool const use_monolithic_scheme);

    //! \name ODESystem interface
    //! @{
    bool isLinear() const override;
    //! @}

    void computeSecondaryVariableConcrete(double const t,
                                          GlobalVector const& x,
                                          int const process_id) override;

private:
    using LocalAssemblerInterface = ThermoHydroMechanicsLocalAssemblerInterface;

    void constructDofTable() override;

    void initializeConcreteProcess(
        NumLib::LocalToGlobalIndexMap const& dof_table,
        MeshLib::Mesh const& mesh,
        unsigned const integration_order) override;

    void assembleConcreteProcess(const double t, GlobalVector const& x,
                                 GlobalMatrix& M, GlobalMatrix& K,
                                 GlobalVector& b) override;

    void assembleWithJacobianConcreteProcess(
        const double t, GlobalVector const& x, GlobalVector const& xdot,
        const double dxdot_dx, const double dx_dx, GlobalMatrix& M,
        GlobalMatrix& K, GlobalVector& b, GlobalMatrix& Jac) override;
    void preTimestepConcreteProcess(GlobalVector const& x, double const t,
                                    double const dt,
                                    const int /*process_id*/) override;

private:
    ThermoHydroMechanicsProcessData<GlobalDim> _process_data;

    std::vector<std::unique_ptr<LocalAssemblerInterface>> _local_assemblers;

    std::vector<MeshLib::Element*> _vec_matrix_elements;
    std::vector<int> _vec_fracture_mat_IDs;
    std::vector<MeshLib::Element*> _vec_all_fracture_elements;
    std::vector<std::vector<MeshLib::Element*>> _vec_fracture_elements;
    std::vector<std::vector<MeshLib::Element*>> _vec_fracture_matrix_elements;
    std::vector<std::vector<MeshLib::Node*>> _vec_fracture_nodes;
    std::vector<std::vector<MeshLib::Node*>> _vec_fracture_nodes_with_tips;
    std::vector<MeshLib::Node*> _vec_junction_nodes;
    std::vector<std::vector<MeshLib::Element*>> _vec_junction_fracture_matrix_elements;

    std::vector<std::unique_ptr<MeshLib::MeshSubset const>> _mesh_subset_fracture_nodes;
    std::unique_ptr<MeshLib::MeshSubset const> _mesh_subset_junction_nodes;
    std::unique_ptr<MeshLib::MeshSubset const> _mesh_subset_matrix_nodes;

    std::vector<MeshLib::Node*> _mesh_nodes_p;
    std::unique_ptr<MeshLib::MeshSubset const> _mesh_subset_nodes_p;
    std::vector<MeshLib::Node*> _mesh_nodes_T;
    std::unique_ptr<MeshLib::MeshSubset const> _mesh_subset_nodes_T;
};

extern template class ThermoHydroMechanicsProcess<2>;
extern template class ThermoHydroMechanicsProcess<3>;

}  // namespace ThermoHydroMechanics
}  // namespace LIE
}  // namespace ProcessLib