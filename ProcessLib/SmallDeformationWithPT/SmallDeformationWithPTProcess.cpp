/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "SmallDeformationWithPTProcess.h"

#include <cassert>

#include "BaseLib/Functional.h"
#include "ProcessLib/SmallDeformation/CreateLocalAssemblers.h"

#include "SmallDeformationWithPTFEM.h"

namespace ProcessLib
{
namespace SmallDeformationWithPT
{
template <int DisplacementDim>
SmallDeformationWithPTProcess<DisplacementDim>::SmallDeformationWithPTProcess(
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    unsigned const integration_order,
    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>&&
        process_variables,
    SmallDeformationWithPTProcessData<DisplacementDim>&& process_data,
    SecondaryVariableCollection&& secondary_variables,
    NumLib::NamedFunctionCaller&& named_function_caller)
    : Process(mesh, std::move(jacobian_assembler), parameters,
              integration_order, std::move(process_variables),
              std::move(secondary_variables), std::move(named_function_caller),
              true),
      _process_data(std::move(process_data))
{
    _integration_point_writer.emplace_back(
        std::make_unique<SigmaIntegrationPointWriter>(
            "sigma_ip",
            static_cast<int>(mesh.getDimension() == 2 ? 4 : 6) /*n components*/,
            2 /*integration order*/, [this]() {
                // Result containing integration point data for each local
                // assembler.
                std::vector<std::vector<double>> result;
                result.resize(_local_assemblers.size());

                for (std::size_t i = 0; i < _local_assemblers.size(); ++i)
                {
                    auto const& local_asm = *_local_assemblers[i];

                    result[i] = local_asm.getSigma();
                }

                return result;
            }));

    _integration_point_writer.emplace_back(
        std::make_unique<SigmaIntegrationPointWriter>(
            "epsilon_ip",
            static_cast<int>(mesh.getDimension() == 2 ? 4 : 6) /*n components*/,
            2 /*integration order*/, [this]() {
                // Result containing integration point data for each local
                // assembler.
                std::vector<std::vector<double>> result;
                result.resize(_local_assemblers.size());

                for (std::size_t i = 0; i < _local_assemblers.size(); ++i)
                {
                    auto const& local_asm = *_local_assemblers[i];

                    result[i] = local_asm.getEpsilon();
                }

                return result;
            }));

    _integration_point_writer.emplace_back(
        std::make_unique<SigmaIntegrationPointWriter>(
            "epsilon_m_ip",
            static_cast<int>(mesh.getDimension() == 2 ? 4 : 6) /*n components*/,
            2 /*integration order*/, [this]() {
                // Result containing integration point data for each local
                // assembler.
                std::vector<std::vector<double>> result;
                result.resize(_local_assemblers.size());

                for (std::size_t i = 0; i < _local_assemblers.size(); ++i)
                {
                    auto const& local_asm = *_local_assemblers[i];

                    result[i] = local_asm.getEpsilonMechanical();
                }

                return result;
            }));
}

template <int DisplacementDim>
bool SmallDeformationWithPTProcess<DisplacementDim>::isLinear() const
{
    return false;
}

template <int DisplacementDim>
void SmallDeformationWithPTProcess<DisplacementDim>::initializeConcreteProcess(
    NumLib::LocalToGlobalIndexMap const& dof_table,
    MeshLib::Mesh const& mesh,
    unsigned const integration_order)
{
    ProcessLib::SmallDeformation::createLocalAssemblers<
        DisplacementDim, SmallDeformationWithPTLocalAssembler>(
        mesh.getElements(), dof_table, _local_assemblers,
        mesh.isAxiallySymmetric(), integration_order, _process_data);

    // TODO move the two data members somewhere else.
    // for extrapolation of secondary variables
    std::vector<MeshLib::MeshSubset> all_mesh_subsets_single_component{
        *_mesh_subset_all_nodes};
    _local_to_global_index_map_single_component.reset(
        new NumLib::LocalToGlobalIndexMap(
            std::move(all_mesh_subsets_single_component),
            // by location order is needed for output
            NumLib::ComponentOrder::BY_LOCATION));

    _secondary_variables.addSecondaryVariable(
        "sigma",
        makeExtrapolator(
            MathLib::KelvinVector::KelvinVectorType<
                DisplacementDim>::RowsAtCompileTime,
            getExtrapolator(), _local_assemblers,
            &SmallDeformationWithPTLocalAssemblerInterface::getIntPtSigma));

    _secondary_variables.addSecondaryVariable(
        "epsilon",
        makeExtrapolator(
            MathLib::KelvinVector::KelvinVectorType<
                DisplacementDim>::RowsAtCompileTime,
            getExtrapolator(), _local_assemblers,
            &SmallDeformationWithPTLocalAssemblerInterface::getIntPtEpsilon));

    auto mesh_prop_pressure_prev = MeshLib::getOrCreateMeshProperty<double>(
        const_cast<MeshLib::Mesh&>(mesh), "pressure_prev",
        MeshLib::MeshItemType::Node, 1);
    mesh_prop_pressure_prev->resize(mesh.getNumberOfNodes());

    auto mesh_prop_pressure = MeshLib::getOrCreateMeshProperty<double>(
        const_cast<MeshLib::Mesh&>(mesh), "pressure",
        MeshLib::MeshItemType::Node, 1);
    mesh_prop_pressure->resize(mesh.getNumberOfNodes());

    auto mesh_prop_temperature_prev = MeshLib::getOrCreateMeshProperty<double>(
        const_cast<MeshLib::Mesh&>(mesh), "temperature_prev",
        MeshLib::MeshItemType::Node, 1);
    mesh_prop_temperature_prev->resize(mesh.getNumberOfNodes());

    auto mesh_prop_temperature = MeshLib::getOrCreateMeshProperty<double>(
        const_cast<MeshLib::Mesh&>(mesh), "temperature",
        MeshLib::MeshItemType::Node, 1);
    mesh_prop_temperature->resize(mesh.getNumberOfNodes());

    for (std::size_t i=0; i<mesh.getNumberOfNodes(); i++)
    {
        SpatialPosition x_pos;
        x_pos.setNodeID(i);
        (*mesh_prop_pressure_prev)[i] = _process_data.p0(0, x_pos)[0];
        (*mesh_prop_pressure)[i] = _process_data.p1(0, x_pos)[0];
        (*mesh_prop_temperature_prev)[i] = _process_data.T0(0, x_pos)[0];
        (*mesh_prop_temperature)[i] = _process_data.T1(0, x_pos)[0];
    }

    // Set initial conditions for integration point data.
    for (auto const& ip_writer : _integration_point_writer)
    {
        // Find the mesh property with integration point writer's name.
        auto const& name = ip_writer->name();
        if (!mesh.getProperties().existsPropertyVector<double>(name))
        {
            continue;
        }
        auto const& mesh_property =
            *mesh.getProperties().template getPropertyVector<double>(name);

        // The mesh property must be defined on integration points.
        if (mesh_property.getMeshItemType() !=
            MeshLib::MeshItemType::IntegrationPoint)
        {
            continue;
        }

        auto const ip_meta_data = getIntegrationPointMetaData(mesh, name);

        // Check the number of components.
        if (ip_meta_data.n_components != mesh_property.getNumberOfComponents())
        {
            OGS_FATAL(
                "Different number of components in meta data (%d) than in "
                "the integration point field data for '%s': %d.",
                ip_meta_data.n_components, name.c_str(),
                mesh_property.getNumberOfComponents());
        }

        // Now we have a properly named vtk's field data array and the
        // corresponding meta data.
        std::size_t position = 0;
        for (auto& local_asm : _local_assemblers)
        {
            std::size_t const integration_points_read =
                local_asm->setIPDataInitialConditions(
                    name, &mesh_property[position],
                    ip_meta_data.integration_order);
            if (integration_points_read == 0)
            {
                OGS_FATAL(
                    "No integration points read in the integration point "
                    "initial conditions set function.");
            }
            position += integration_points_read * ip_meta_data.n_components;
        }
    }
}

template <int DisplacementDim>
void SmallDeformationWithPTProcess<DisplacementDim>::assembleConcreteProcess(
    const double t, GlobalVector const& x, GlobalMatrix& M, GlobalMatrix& K,
    GlobalVector& b)
{
    DBUG("Assemble SmallDeformationWithPTProcess.");

    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
       dof_table = {std::ref(*_local_to_global_index_map)};
    const int process_id = 0;
    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];

    // Call global assembler for each local assembly item.
    GlobalExecutor::executeSelectedMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assemble, _local_assemblers,
        pv.getActiveElementIDs(), dof_table, t, x, M, K, b,
        _coupled_solutions);
}

template <int DisplacementDim>
void SmallDeformationWithPTProcess<DisplacementDim>::
    assembleWithJacobianConcreteProcess(const double t, GlobalVector const& x,
                                        GlobalVector const& xdot,
                                        const double dxdot_dx,
                                        const double dx_dx, GlobalMatrix& M,
                                        GlobalMatrix& K, GlobalVector& b,
                                        GlobalMatrix& Jac)
{
    DBUG("AssembleJacobian SmallDeformationWithPTProcess.");

    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
       dof_table = {std::ref(*_local_to_global_index_map)};
     const int process_id = 0;
    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];

    // Call global assembler for each local assembly item.
    GlobalExecutor::executeSelectedMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assembleWithJacobian,
        _local_assemblers, pv.getActiveElementIDs(), dof_table, t, x,
        xdot, dxdot_dx, dx_dx, M, K, b, Jac, _coupled_solutions);
}

template <int DisplacementDim>
void SmallDeformationWithPTProcess<DisplacementDim>::preTimestepConcreteProcess(
    GlobalVector const& x, double const t, double const dt,
    const int process_id)
{
    DBUG("PreTimestep SmallDeformationWithPTProcess.");

    _process_data.dt = dt;
    _process_data.t = t;

    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];

    GlobalExecutor::executeSelectedMemberOnDereferenced(
        &SmallDeformationWithPTLocalAssemblerInterface::preTimestep, _local_assemblers,
        pv.getActiveElementIDs(), *_local_to_global_index_map, x, t,
        dt);
}

template <int DisplacementDim>
void SmallDeformationWithPTProcess<DisplacementDim>::postTimestepConcreteProcess(
    GlobalVector const& x, const double /*t*/, const double /*delta_t*/,
    int const process_id)
{
    DBUG("PostTimestep SmallDeformationWithPTProcess.");

    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];

    GlobalExecutor::executeSelectedMemberOnDereferenced(
        &SmallDeformationWithPTLocalAssemblerInterface::postTimestep,
        _local_assemblers, pv.getActiveElementIDs(),
        *_local_to_global_index_map, x);
}

template class SmallDeformationWithPTProcess<2>;
template class SmallDeformationWithPTProcess<3>;

}  // namespace SmallDeformationWithPT
}  // namespace ProcessLib
