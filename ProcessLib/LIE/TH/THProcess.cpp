
#include "THProcess.h"

#include "MeshLib/ElementCoordinatesMappingLocal.h"
#include "MeshLib/ElementStatus.h"
#include "MeshLib/Elements/Utils.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/MeshInformation.h"
#include "MeshLib/Properties.h"

#include "NumLib/DOF/DOFTableUtil.h"
#include "NumLib/DOF/LocalToGlobalIndexMap.h"

#include "ParameterLib/MeshElementParameter.h"

#include "ProcessLib/LIE/Common/BranchProperty.h"
#include "ProcessLib/LIE/Common/JunctionProperty.h"
#include "ProcessLib/LIE/Common/MeshUtils.h"
#include "ProcessLib/LIE/TH/LocalAssembler/CreateLocalAssemblers.h"
#include "ProcessLib/LIE/TH/LocalAssembler/THLocalAssemblerFracture.h"
#include "ProcessLib/LIE/TH/LocalAssembler/THLocalAssemblerMatrix.h"

namespace ProcessLib
{
namespace LIE
{
namespace TH
{
template <int GlobalDim>
THProcess<GlobalDim>::THProcess(
    std::string name,
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    unsigned const integration_order,
    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>&&
        process_variables,
    THProcessData<GlobalDim>&& process_data,
    SecondaryVariableCollection&& secondary_variables,
    NumLib::NamedFunctionCaller&& named_function_caller,
    bool const use_monolithic_scheme)
    : Process(std::move(name), mesh, std::move(jacobian_assembler), parameters,
              integration_order, std::move(process_variables),
              std::move(secondary_variables), std::move(named_function_caller),
              use_monolithic_scheme),
      _process_data(std::move(process_data))
{
    INFO("[LIE/TH] looking for fracture elements in the given mesh");
    std::vector<std::pair<std::size_t, std::vector<int>>>
        vec_branch_nodeID_matIDs;
    std::vector<std::pair<std::size_t, std::vector<int>>>
        vec_junction_nodeID_matIDs;
    getFractureMatrixDataInMesh(mesh, _vec_matrix_elements,
                                _vec_fracture_mat_IDs, _vec_fracture_elements,
                                _vec_fracture_matrix_elements,
                                _vec_fracture_nodes, _vec_fracture_nodes_with_tips, vec_branch_nodeID_matIDs,
                                vec_junction_nodeID_matIDs);

    if (_vec_fracture_mat_IDs.size() !=
        _process_data.fracture_properties.size())
    {
        OGS_FATAL(
            "The number of the given fracture properties (%d) are not "
            "consistent"
            " with the number of fracture groups in a mesh (%d).",
            _process_data.fracture_properties.size(),
            _vec_fracture_mat_IDs.size());
    }

    // create a map from a material ID to a fracture ID
    auto max_frac_mat_id = std::max_element(_vec_fracture_mat_IDs.begin(),
                                            _vec_fracture_mat_IDs.end());
    _process_data.map_materialID_to_fractureID.resize(*max_frac_mat_id + 1);
    for (unsigned i = 0; i < _vec_fracture_mat_IDs.size(); i++)
    {
        _process_data.map_materialID_to_fractureID[_vec_fracture_mat_IDs[i]] =
            i;
    }

    // create a table of connected fracture IDs for each element
    _process_data.vec_ele_connected_fractureIDs.resize(
        mesh.getNumberOfElements());
    for (unsigned i = 0; i < _vec_fracture_matrix_elements.size(); i++)
    {
        for (auto e : _vec_fracture_matrix_elements[i])
        {
            _process_data.vec_ele_connected_fractureIDs[e->getID()].push_back(
                i);
        }
    }

    // set fracture property
    for (auto& fracture_prop : _process_data.fracture_properties)
    {
        // based on the 1st element assuming a fracture forms a straight line
        setFractureProperty(
            GlobalDim,
            *_vec_fracture_elements[fracture_prop->fracture_id][0],
            *fracture_prop);
    }

    // set branches
    for (auto& vec_branch_nodeID_matID : vec_branch_nodeID_matIDs)
    {
        auto master_matId = vec_branch_nodeID_matID.second[0];
        auto slave_matId = vec_branch_nodeID_matID.second[1];
        auto& master_frac =
            *_process_data.fracture_properties
                [_process_data.map_materialID_to_fractureID[master_matId]];
        auto& slave_frac =
            *_process_data.fracture_properties
                [_process_data.map_materialID_to_fractureID[slave_matId]];

        master_frac.branches_master.push_back(
            createBranchProperty(*mesh.getNode(vec_branch_nodeID_matID.first),
                                 master_frac, slave_frac));

        slave_frac.branches_slave.push_back(
            createBranchProperty(*mesh.getNode(vec_branch_nodeID_matID.first),
                                 master_frac, slave_frac));
    }

    // set junctions
    for (auto& vec_junction_nodeID_matID : vec_junction_nodeID_matIDs)
    {
        _vec_junction_nodes.push_back(const_cast<MeshLib::Node*>(
            _mesh.getNode(vec_junction_nodeID_matID.first)));
    }
    for (std::size_t i = 0; i < vec_junction_nodeID_matIDs.size(); i++)
    {
        auto const& material_ids = vec_junction_nodeID_matIDs[i].second;
        assert(material_ids.size() == 2);
        std::array<int, 2> fracture_ids{
            {_process_data.map_materialID_to_fractureID[material_ids[0]],
             _process_data.map_materialID_to_fractureID[material_ids[1]]}};

        _process_data.junction_properties.emplace_back(
            i, *mesh.getNode(vec_junction_nodeID_matIDs[i].first),
            fracture_ids);
    }

    // create a table of connected junction IDs for each element
    _process_data.vec_ele_connected_junctionIDs.resize(
        mesh.getNumberOfElements());
    for (unsigned i = 0; i < vec_junction_nodeID_matIDs.size(); i++)
    {
        auto node = mesh.getNode(vec_junction_nodeID_matIDs[i].first);
        for (auto e : node->getElements())
        {
            _process_data.vec_ele_connected_junctionIDs[e->getID()].push_back(
                i);
        }
    }

    // create a table of junction node and connected elements
    _vec_junction_fracture_matrix_elements.resize(
        vec_junction_nodeID_matIDs.size());
    for (unsigned i = 0; i < vec_junction_nodeID_matIDs.size(); i++)
    {
        auto node = mesh.getNode(vec_junction_nodeID_matIDs[i].first);
        for (auto e : node->getElements())
        {
            _vec_junction_fracture_matrix_elements[i].push_back(e);
        }
    }

    //
    // If Neumann BCs for the displacement_jump variable are required they need
    // special treatment because of the levelset function. The implementation
    // exists in the version 6.1.0 (e54815cc07ee89c81f953a4955b1c788595dd725)
    // and was removed due to lack of applications.
    //

    if (!_process_data.deactivate_matrix_in_flow)
    {
        _process_data.p_element_status =
            std::make_unique<MeshLib::ElementStatus>(&mesh);
    }
    else
    {
        auto const range =
            MeshLib::MeshInformation::getValueBounds<int>(mesh, "MaterialIDs");
        if (!range)
        {
            OGS_FATAL(
                "Could not get minimum/maximum ranges values for the "
                "MaterialIDs property in the mesh '%s'.",
                mesh.getName().c_str());
        }

        std::vector<int> vec_p_inactive_matIDs;
        for (int matID = range->first; matID <= range->second; matID++)
        {
            if (std::find(_vec_fracture_mat_IDs.begin(),
                          _vec_fracture_mat_IDs.end(),
                          matID) == _vec_fracture_mat_IDs.end())
            {
                vec_p_inactive_matIDs.push_back(matID);
            }
        }
        _process_data.p_element_status =
            std::make_unique<MeshLib::ElementStatus>(&mesh,
                                                     vec_p_inactive_matIDs);

        const int monolithic_process_id = 0;
        ProcessVariable const& pv_p =
            getProcessVariables(monolithic_process_id)[0];
        _process_data.p0 = &pv_p.getInitialCondition();
    }

    MeshLib::PropertyVector<int> const* material_ids(
        mesh.getProperties().getPropertyVector<int>("MaterialIDs"));
    _process_data.mesh_prop_materialIDs = material_ids;
}

template <int GlobalDim>
void THProcess<GlobalDim>::constructDofTable()
{
    //------------------------------------------------------------
    // prepare mesh subsets to define DoFs
    //------------------------------------------------------------
    // for extrapolation
    _mesh_subset_all_nodes =
        std::make_unique<MeshLib::MeshSubset>(_mesh, _mesh.getNodes());
    // pressure
    _mesh_nodes_p = MeshLib::getBaseNodes(
        _process_data.p_element_status->getActiveElements());
    _mesh_subset_nodes_p =
        std::make_unique<MeshLib::MeshSubset>(_mesh, _mesh_nodes_p);
    // temperature
    _mesh_nodes_T = MeshLib::getBaseNodes(_mesh.getElements());
    _mesh_subset_nodes_T =
        std::make_unique<MeshLib::MeshSubset>(_mesh, _mesh_nodes_T);

    // Collect the mesh subsets in a vector.
    std::vector<MeshLib::MeshSubset> all_mesh_subsets;
    std::vector<int> vec_n_components;
    std::vector<std::vector<MeshLib::Element*> const*> vec_var_elements;
    // pressure
    vec_n_components.push_back(1);
    all_mesh_subsets.emplace_back(*_mesh_subset_nodes_p);
    if (!_process_data.deactivate_matrix_in_flow)
    {
        vec_var_elements.push_back(&_mesh.getElements());
    }
    else
    {
        // TODO set elements including active nodes for pressure.
        // cannot use ElementStatus
        for (auto& vec : _vec_fracture_matrix_elements)
            vec_var_elements.push_back(&vec);
    }
    // temperature
    vec_n_components.push_back(1);
    all_mesh_subsets.emplace_back(*_mesh_subset_nodes_T);
    vec_var_elements.push_back(&_mesh.getElements());

    INFO("[LIE/TH] creating a DoF table");
    _local_to_global_index_map =
        std::make_unique<NumLib::LocalToGlobalIndexMap>(
            std::move(all_mesh_subsets),
            vec_n_components,
            vec_var_elements,
            NumLib::ComponentOrder::BY_COMPONENT);

    DBUG("created %d DoF", _local_to_global_index_map->size());
}

template <int GlobalDim>
void THProcess<GlobalDim>::initializeConcreteProcess(
    NumLib::LocalToGlobalIndexMap const& dof_table,
    MeshLib::Mesh const& mesh,
    unsigned const integration_order)
{
    assert(mesh.getDimension() == GlobalDim);
    INFO("[LIE/TH] creating local assemblers");
    const int monolithic_process_id = 0;
    ProcessLib::LIE::TH::createLocalAssemblers<
        GlobalDim, THLocalAssemblerMatrix,
        THLocalAssemblerFracture>(
        mesh.getElements(), dof_table,
        // use displacment process variable for shapefunction order
        getProcessVariables(
            monolithic_process_id)[0].get().getShapeFunctionOrder(),
            _local_assemblers, mesh.isAxiallySymmetric(), integration_order,
            _process_data);

    auto mesh_prop_velocity = MeshLib::getOrCreateMeshProperty<double>(
        const_cast<MeshLib::Mesh&>(mesh), "velocity",
        MeshLib::MeshItemType::Cell, 3);
    mesh_prop_velocity->resize(mesh.getNumberOfElements() * 3);
    _process_data.mesh_prop_velocity = mesh_prop_velocity;

    _process_data.mesh_prop_hydraulic_flow =
        MeshLib::getOrCreateMeshProperty<double>(
            const_cast<MeshLib::Mesh&>(mesh), "HydraulicFlow",
            MeshLib::MeshItemType::Node, 1);
    assert(_process_data.mesh_prop_hydraulic_flow->size() ==
            mesh.getNumberOfNodes());

    _process_data.mesh_prop_thermal_flow =
        MeshLib::getOrCreateMeshProperty<double>(
            const_cast<MeshLib::Mesh&>(mesh), "ThermalFlow",
            MeshLib::MeshItemType::Node, 1);
    assert(_process_data.mesh_prop_thermal_flow->size() ==
            mesh.getNumberOfNodes());
}

template <int GlobalDim>
void THProcess<GlobalDim>::computeSecondaryVariableConcrete(
    const double t, GlobalVector const& x, int const process_id)
{
    DBUG("Compute the secondary variables for THProcess.");
    const auto& dof_table = getDOFTable(process_id);

    {
        ProcessLib::ProcessVariable const& pv =
            getProcessVariables(process_id)[0];

        GlobalExecutor::executeSelectedMemberOnDereferenced(
            &THLocalAssemblerInterface::computeSecondaryVariable,
            _local_assemblers, pv.getActiveElementIDs(),
            dof_table, t, x, _coupled_solutions);
    }


    MathLib::LinAlg::setLocalAccessibleVector(x);
}

template <int GlobalDim>
bool THProcess<GlobalDim>::isLinear() const
{
    return false;
}

template <int GlobalDim>
void THProcess<GlobalDim>::assembleConcreteProcess(
    const double t, GlobalVector const& x, GlobalMatrix& M, GlobalMatrix& K,
    GlobalVector& b)
{
    DBUG("Assemble THProcess.");

    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
        dof_table = {std::ref(*_local_to_global_index_map)};
    // Call global assembler for each local assembly item.
    GlobalExecutor::executeMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assemble, _local_assemblers,
        dof_table, t, x, M, K, b, _coupled_solutions);
}

template <int GlobalDim>
void THProcess<GlobalDim>::assembleWithJacobianConcreteProcess(
    const double t, GlobalVector const& x, GlobalVector const& xdot,
    const double dxdot_dx, const double dx_dx, GlobalMatrix& M, GlobalMatrix& K,
    GlobalVector& b, GlobalMatrix& Jac)
{
    DBUG("AssembleWithJacobian THProcess.");

    const int process_id =
        _use_monolithic_scheme ? 0 : _coupled_solutions->process_id;
    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];

    // Call global assembler for each local assembly item.
    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
       dof_table = {std::ref(*_local_to_global_index_map)};
    GlobalExecutor::executeSelectedMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assembleWithJacobian,
        _local_assemblers, pv.getActiveElementIDs(), dof_table, t, x,
        xdot, dxdot_dx, dx_dx, M, K, b, Jac, _coupled_solutions);

    // b.write("b.txt");
    // Jac.write("J.txt");

    auto copyRhs = [&](int const variable_id, auto& output_vector) {
        transformVariableFromGlobalVector(b, variable_id,
                                          *_local_to_global_index_map,
                                          output_vector, std::negate<double>());
    };
    copyRhs(0, *_process_data.mesh_prop_hydraulic_flow);
    copyRhs(1, *_process_data.mesh_prop_thermal_flow);
}

template <int GlobalDim>
void THProcess<GlobalDim>::preTimestepConcreteProcess(
    GlobalVector const& x, double const t, double const dt,
    const int process_id)
{
    DBUG("PreTimestep THProcess.");

    _process_data.dt = dt;
    _process_data.t = t;

    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];

    GlobalExecutor::executeSelectedMemberOnDereferenced(
        &THLocalAssemblerInterface::preTimestep, _local_assemblers,
        pv.getActiveElementIDs(), *_local_to_global_index_map,
        x, t, dt);
}

// ------------------------------------------------------------------------------------
// template instantiation
// ------------------------------------------------------------------------------------
template class THProcess<1>;
template class THProcess<2>;
template class THProcess<3>;

}  // namespace TH
}  // namespace LIE
}  // namespace ProcessLib
