
#include "ThermoHydroMechanicsProcess.h"

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
#include "ProcessLib/LIE/ThermoHydroMechanics/LocalAssembler/CreateLocalAssemblers.h"
#include "ProcessLib/LIE/ThermoHydroMechanics/LocalAssembler/ThermoHydroMechanicsLocalAssemblerFracture.h"
#include "ProcessLib/LIE/ThermoHydroMechanics/LocalAssembler/ThermoHydroMechanicsLocalAssemblerMatrix.h"
#include "ProcessLib/LIE/ThermoHydroMechanics/LocalAssembler/ThermoHydroMechanicsLocalAssemblerMatrixNearFracture.h"

namespace ProcessLib
{
namespace LIE
{
namespace ThermoHydroMechanics
{
template <int GlobalDim>
ThermoHydroMechanicsProcess<GlobalDim>::ThermoHydroMechanicsProcess(
    std::string name,
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    unsigned const integration_order,
    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>&&
        process_variables,
    ThermoHydroMechanicsProcessData<GlobalDim>&& process_data,
    SecondaryVariableCollection&& secondary_variables,
    NumLib::NamedFunctionCaller&& named_function_caller,
    bool const use_monolithic_scheme)
    : Process(std::move(name), mesh, std::move(jacobian_assembler), parameters,
              integration_order, std::move(process_variables),
              std::move(secondary_variables), std::move(named_function_caller),
              use_monolithic_scheme),
      _process_data(std::move(process_data))
{
    INFO("[LIE/THM] looking for fracture elements in the given mesh");
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
    // need to use a custom Neumann BC assembler for displacement jumps
    for (ProcessVariable& pv : getProcessVariables())
    {
        if (pv.getName().find("displacement_jump") == std::string::npos)
            continue;
        pv.setBoundaryConditionBuilder(
                    std::unique_ptr<ProcessLib::BoundaryConditionBuilder>(
                        new BoundaryConditionBuilder(*_process_data.fracture_property.get())));
    }

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
void ThermoHydroMechanicsProcess<GlobalDim>::constructDofTable()
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
    // regular u
    _mesh_subset_matrix_nodes =
        std::make_unique<MeshLib::MeshSubset>(_mesh, _mesh.getNodes());
    // u jump
    for (unsigned i = 0; i < _vec_fracture_nodes.size(); i++)
    {
        _mesh_subset_fracture_nodes.push_back(
            std::make_unique<MeshLib::MeshSubset const>(
                _mesh, _vec_fracture_nodes[i]));
    }
    // enrichment for junctions
    _mesh_subset_junction_nodes =
        std::make_unique<MeshLib::MeshSubset>(_mesh, _vec_junction_nodes);

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
        for (auto& vec : _vec_fracture_matrix_elements)
            std::copy(vec.begin(), vec.end(),
                    std::back_inserter(_vec_all_fracture_elements));

        BaseLib::makeVectorUnique(_vec_all_fracture_elements);
        vec_var_elements.push_back(&_vec_all_fracture_elements);
    }
    // temperature
    vec_n_components.push_back(1);
    all_mesh_subsets.emplace_back(*_mesh_subset_nodes_T);
    vec_var_elements.push_back(&_mesh.getElements());
    // regular displacement
    vec_n_components.push_back(GlobalDim);
    std::generate_n(std::back_inserter(all_mesh_subsets), GlobalDim,
                    [&]() { return *_mesh_subset_matrix_nodes; });
    for (auto& ms : _mesh_subset_fracture_nodes)
    {
        vec_n_components.push_back(GlobalDim);
        std::generate_n(std::back_inserter(all_mesh_subsets),
                        GlobalDim,
                        [&]() { return *ms; });
    }
    if (!_vec_junction_nodes.empty())
    {
        std::generate_n(std::back_inserter(all_mesh_subsets),
                        GlobalDim,
                        [&]() { return *_mesh_subset_junction_nodes; });
        for (unsigned i=0; i<_vec_junction_nodes.size(); i++)
        vec_n_components.push_back(GlobalDim);
    }

    vec_var_elements.push_back(&_vec_matrix_elements);
    for (unsigned i = 0; i < _vec_fracture_matrix_elements.size(); i++)
    {
        vec_var_elements.push_back(&_vec_fracture_matrix_elements[i]);
    }
    for (unsigned i = 0; i < _vec_junction_fracture_matrix_elements.size(); i++)
    {
        vec_var_elements.push_back(&_vec_junction_fracture_matrix_elements[i]);
    }

    INFO("[LIE/THM] creating a DoF table");
    _local_to_global_index_map =
        std::make_unique<NumLib::LocalToGlobalIndexMap>(
            std::move(all_mesh_subsets),
            vec_n_components,
            vec_var_elements,
            NumLib::ComponentOrder::BY_COMPONENT);

    DBUG("created %d DoF", _local_to_global_index_map->size());
}

template <int GlobalDim>
void ThermoHydroMechanicsProcess<GlobalDim>::initializeConcreteProcess(
    NumLib::LocalToGlobalIndexMap const& dof_table,
    MeshLib::Mesh const& mesh,
    unsigned const integration_order)
{
    assert(mesh.getDimension() == GlobalDim);
    INFO("[LIE/THM] creating local assemblers");
    const int monolithic_process_id = 0;
    const int monolithic_pv_id_u = 2;
    ProcessLib::LIE::ThermoHydroMechanics::createLocalAssemblers<
        GlobalDim, ThermoHydroMechanicsLocalAssemblerMatrix,
        ThermoHydroMechanicsLocalAssemblerMatrixNearFracture,
        ThermoHydroMechanicsLocalAssemblerFracture>(
        mesh.getElements(), dof_table,
        // use displacment process variable for shapefunction order
        getProcessVariables(
            monolithic_process_id)[monolithic_pv_id_u].get().getShapeFunctionOrder(),
            _local_assemblers, mesh.isAxiallySymmetric(), integration_order,
            _process_data);

    auto mesh_prop_sigma_xx = MeshLib::getOrCreateMeshProperty<double>(
        const_cast<MeshLib::Mesh&>(mesh), "stress_xx",
        MeshLib::MeshItemType::Cell, 1);
    mesh_prop_sigma_xx->resize(mesh.getNumberOfElements());
    _process_data.mesh_prop_stress_xx = mesh_prop_sigma_xx;

    auto mesh_prop_sigma_yy = MeshLib::getOrCreateMeshProperty<double>(
        const_cast<MeshLib::Mesh&>(mesh), "stress_yy",
        MeshLib::MeshItemType::Cell, 1);
    mesh_prop_sigma_yy->resize(mesh.getNumberOfElements());
    _process_data.mesh_prop_stress_yy = mesh_prop_sigma_yy;

    auto mesh_prop_sigma_zz = MeshLib::getOrCreateMeshProperty<double>(
        const_cast<MeshLib::Mesh&>(mesh), "stress_zz",
        MeshLib::MeshItemType::Cell, 1);
    mesh_prop_sigma_zz->resize(mesh.getNumberOfElements());
    _process_data.mesh_prop_stress_zz = mesh_prop_sigma_zz;

    auto mesh_prop_sigma_xy = MeshLib::getOrCreateMeshProperty<double>(
        const_cast<MeshLib::Mesh&>(mesh), "stress_xy",
        MeshLib::MeshItemType::Cell, 1);
    mesh_prop_sigma_xy->resize(mesh.getNumberOfElements());
    _process_data.mesh_prop_stress_xy = mesh_prop_sigma_xy;

    if (GlobalDim == 3)
    {
        auto mesh_prop_sigma_xz = MeshLib::getOrCreateMeshProperty<double>(
            const_cast<MeshLib::Mesh&>(mesh), "stress_xz",
            MeshLib::MeshItemType::Cell, 1);
        mesh_prop_sigma_xz->resize(mesh.getNumberOfElements());
        _process_data.mesh_prop_stress_xz = mesh_prop_sigma_xz;

        auto mesh_prop_sigma_yz = MeshLib::getOrCreateMeshProperty<double>(
            const_cast<MeshLib::Mesh&>(mesh), "stress_yz",
            MeshLib::MeshItemType::Cell, 1);
        mesh_prop_sigma_yz->resize(mesh.getNumberOfElements());
        _process_data.mesh_prop_stress_yz = mesh_prop_sigma_yz;
    }

    auto mesh_prop_epsilon_xx = MeshLib::getOrCreateMeshProperty<double>(
        const_cast<MeshLib::Mesh&>(mesh), "strain_xx",
        MeshLib::MeshItemType::Cell, 1);
    mesh_prop_epsilon_xx->resize(mesh.getNumberOfElements());
    _process_data.mesh_prop_strain_xx = mesh_prop_epsilon_xx;

    auto mesh_prop_epsilon_yy = MeshLib::getOrCreateMeshProperty<double>(
        const_cast<MeshLib::Mesh&>(mesh), "strain_yy",
        MeshLib::MeshItemType::Cell, 1);
    mesh_prop_epsilon_yy->resize(mesh.getNumberOfElements());
    _process_data.mesh_prop_strain_yy = mesh_prop_epsilon_yy;

    auto mesh_prop_epsilon_zz = MeshLib::getOrCreateMeshProperty<double>(
        const_cast<MeshLib::Mesh&>(mesh), "strain_zz",
        MeshLib::MeshItemType::Cell, 1);
    mesh_prop_epsilon_zz->resize(mesh.getNumberOfElements());
    _process_data.mesh_prop_strain_zz = mesh_prop_epsilon_zz;

    auto mesh_prop_epsilon_xy = MeshLib::getOrCreateMeshProperty<double>(
        const_cast<MeshLib::Mesh&>(mesh), "strain_xy",
        MeshLib::MeshItemType::Cell, 1);
    mesh_prop_epsilon_xy->resize(mesh.getNumberOfElements());
    _process_data.mesh_prop_strain_xy = mesh_prop_epsilon_xy;

    if (GlobalDim == 3)
    {
        auto mesh_prop_epsilon_xz = MeshLib::getOrCreateMeshProperty<double>(
            const_cast<MeshLib::Mesh&>(mesh), "strain_xz",
            MeshLib::MeshItemType::Cell, 1);
        mesh_prop_epsilon_xz->resize(mesh.getNumberOfElements());
        _process_data.mesh_prop_strain_xz = mesh_prop_epsilon_xz;

        auto mesh_prop_epsilon_yz = MeshLib::getOrCreateMeshProperty<double>(
            const_cast<MeshLib::Mesh&>(mesh), "strain_yz",
            MeshLib::MeshItemType::Cell, 1);
        mesh_prop_epsilon_yz->resize(mesh.getNumberOfElements());
        _process_data.mesh_prop_strain_yz = mesh_prop_epsilon_yz;
    }

    auto mesh_prop_velocity = MeshLib::getOrCreateMeshProperty<double>(
        const_cast<MeshLib::Mesh&>(mesh), "velocity",
        MeshLib::MeshItemType::Cell, 3);
    mesh_prop_velocity->resize(mesh.getNumberOfElements() * 3);
    _process_data.mesh_prop_velocity = mesh_prop_velocity;

    auto const n_enrich_var = _process_data.fracture_properties.size() + _process_data.junction_properties.size();

    if (!_vec_fracture_elements.empty())
    {
        for (MeshLib::Element const* e : _mesh.getElements())
        {
            if (e->getDimension() < GlobalDim)
            {
                continue;
            }

            Eigen::Vector3d const pt(e->getCenterOfGravity().getCoords());
            std::vector<FractureProperty*> e_fracture_props;
            std::unordered_map<int, int> e_fracID_to_local;
            unsigned tmpi = 0;
            for (auto fid :
                _process_data.vec_ele_connected_fractureIDs[e->getID()])
            {
                e_fracture_props.push_back(&*_process_data.fracture_properties[fid]);
                e_fracID_to_local.insert({fid, tmpi++});
            }
            std::vector<JunctionProperty*> e_junction_props;
            std::unordered_map<int, int> e_juncID_to_local;
            tmpi = 0;
            for (auto fid :
                _process_data.vec_ele_connected_junctionIDs[e->getID()])
            {
                e_junction_props.push_back(&_process_data.junction_properties[fid]);
                e_juncID_to_local.insert({fid, tmpi++});
            }
            std::vector<double> const levelsets(uGlobalEnrichments(
                e_fracture_props, e_junction_props, e_fracID_to_local, pt));

            for (unsigned i = 0; i < e_fracture_props.size(); i++)
            {
                auto mesh_prop_levelset = MeshLib::getOrCreateMeshProperty<double>(
                    const_cast<MeshLib::Mesh&>(mesh),
                    "levelset" +
                        std::to_string(e_fracture_props[i]->fracture_id + 1),
                    MeshLib::MeshItemType::Cell, 1);
                mesh_prop_levelset->resize(mesh.getNumberOfElements());
                (*mesh_prop_levelset)[e->getID()] = levelsets[i];
            }
            for (unsigned i = 0; i < e_junction_props.size(); i++)
            {
                auto mesh_prop_levelset = MeshLib::getOrCreateMeshProperty<double>(
                    const_cast<MeshLib::Mesh&>(mesh),
                    "levelset" +
                        std::to_string(e_junction_props[i]->junction_id + 1 +
                                    _process_data.fracture_properties.size()),
                    MeshLib::MeshItemType::Cell, 1);
                mesh_prop_levelset->resize(mesh.getNumberOfElements());
                (*mesh_prop_levelset)[e->getID()] =
                    levelsets[i + e_fracture_props.size()];
            }
        }

        auto mesh_prop_w_n = MeshLib::getOrCreateMeshProperty<double>(
            const_cast<MeshLib::Mesh&>(mesh), "w_n",
            MeshLib::MeshItemType::Cell, 1);
        mesh_prop_w_n->resize(mesh.getNumberOfElements());
        auto mesh_prop_w_s = MeshLib::getOrCreateMeshProperty<double>(
            const_cast<MeshLib::Mesh&>(mesh), "w_s",
            MeshLib::MeshItemType::Cell, 1);
        mesh_prop_w_s->resize(mesh.getNumberOfElements());
        _process_data.mesh_prop_w_n = mesh_prop_w_n;
        _process_data.mesh_prop_w_s = mesh_prop_w_s;

        auto mesh_prop_b = MeshLib::getOrCreateMeshProperty<double>(
            const_cast<MeshLib::Mesh&>(mesh), "aperture",
            MeshLib::MeshItemType::Cell, 1);
        mesh_prop_b->resize(mesh.getNumberOfElements());
        auto const& mesh_prop_matid = *_process_data.mesh_prop_materialIDs;
        for (auto const& frac : _process_data.fracture_properties)
        {
            for (MeshLib::Element const* e : _mesh.getElements())
            {
                if (e->getDimension() != GlobalDim-1)
                    continue;
                if (mesh_prop_matid[e->getID()] != frac->mat_id)
                    continue;

                ParameterLib::SpatialPosition x;
                x.setElementID(e->getID());
                (*mesh_prop_b)[e->getID()] = frac->aperture0(0, x)[0];
            }
        }
        _process_data.mesh_prop_b = mesh_prop_b;

        auto mesh_prop_k_f = MeshLib::getOrCreateMeshProperty<double>(
            const_cast<MeshLib::Mesh&>(mesh), "k_f",
            MeshLib::MeshItemType::Cell, 1);
        mesh_prop_k_f->resize(mesh.getNumberOfElements());
        _process_data.mesh_prop_k_f = mesh_prop_k_f;

        auto mesh_prop_fracture_stress_shear =
            MeshLib::getOrCreateMeshProperty<double>(
                const_cast<MeshLib::Mesh&>(mesh), "f_stress_s",
                MeshLib::MeshItemType::Cell, 1);
        mesh_prop_fracture_stress_shear->resize(mesh.getNumberOfElements());
        _process_data.mesh_prop_fracture_stress_shear =
            mesh_prop_fracture_stress_shear;

        auto mesh_prop_fracture_stress_normal =
            MeshLib::getOrCreateMeshProperty<double>(
                const_cast<MeshLib::Mesh&>(mesh), "f_stress_n",
                MeshLib::MeshItemType::Cell, 1);
        mesh_prop_fracture_stress_normal->resize(mesh.getNumberOfElements());
        _process_data.mesh_prop_fracture_stress_normal =
            mesh_prop_fracture_stress_normal;

        auto mesh_prop_fracture_shear_failure =
            MeshLib::getOrCreateMeshProperty<double>(
                const_cast<MeshLib::Mesh&>(mesh), "f_shear_failure",
                MeshLib::MeshItemType::Cell, 1);
        mesh_prop_fracture_shear_failure->resize(mesh.getNumberOfElements());
        _process_data.mesh_prop_fracture_shear_failure =
            mesh_prop_fracture_shear_failure;

        auto mesh_prop_nodal_w = MeshLib::getOrCreateMeshProperty<double>(
            const_cast<MeshLib::Mesh&>(mesh), "nodal_w",
            MeshLib::MeshItemType::Node, GlobalDim);
        mesh_prop_nodal_w->resize(mesh.getNumberOfNodes() * GlobalDim);
        _process_data.mesh_prop_nodal_w = mesh_prop_nodal_w;

        auto mesh_prop_nodal_b = MeshLib::getOrCreateMeshProperty<double>(
            const_cast<MeshLib::Mesh&>(mesh), "nodal_aperture",
            MeshLib::MeshItemType::Node, 1);
        mesh_prop_nodal_b->resize(mesh.getNumberOfNodes());
        _process_data.mesh_prop_nodal_b = mesh_prop_nodal_b;

        if (GlobalDim == 3)
        {
            auto mesh_prop_w_s2 = MeshLib::getOrCreateMeshProperty<double>(
                const_cast<MeshLib::Mesh&>(mesh), "w_s2",
                MeshLib::MeshItemType::Cell, 1);
            mesh_prop_w_s2->resize(mesh.getNumberOfElements());
            _process_data.mesh_prop_w_s2 = mesh_prop_w_s2;

            auto mesh_prop_fracture_stress_shear2 =
                MeshLib::getOrCreateMeshProperty<double>(
                    const_cast<MeshLib::Mesh&>(mesh), "f_stress_s2",
                    MeshLib::MeshItemType::Cell, 1);
            mesh_prop_fracture_stress_shear2->resize(
                mesh.getNumberOfElements());
            _process_data.mesh_prop_fracture_stress_shear2 =
                mesh_prop_fracture_stress_shear2;
        }

        auto mesh_prop_nodal_p = MeshLib::getOrCreateMeshProperty<double>(
            const_cast<MeshLib::Mesh&>(mesh), "pressure_interpolated",
            MeshLib::MeshItemType::Node, 1);
        mesh_prop_nodal_p->resize(mesh.getNumberOfNodes());
        _process_data.mesh_prop_nodal_p = mesh_prop_nodal_p;

        auto mesh_prop_nodal_T = MeshLib::getOrCreateMeshProperty<double>(
            const_cast<MeshLib::Mesh&>(mesh), "temperature_interpolated",
            MeshLib::MeshItemType::Node, 1);
        mesh_prop_nodal_T->resize(mesh.getNumberOfNodes());
        _process_data.mesh_prop_nodal_T = mesh_prop_nodal_T;

        _process_data.mesh_prop_nodal_forces =
            MeshLib::getOrCreateMeshProperty<double>(
                const_cast<MeshLib::Mesh&>(mesh), "NodalForces",
                MeshLib::MeshItemType::Node, GlobalDim);
        assert(_process_data.mesh_prop_nodal_forces->size() ==
               GlobalDim * mesh.getNumberOfNodes());

        _process_data.vec_mesh_prop_nodal_forces_jump.resize(n_enrich_var);
        for (unsigned ig=0; ig<n_enrich_var; ig++)
        {
            _process_data.vec_mesh_prop_nodal_forces_jump[ig] =
                MeshLib::getOrCreateMeshProperty<double>(
                    const_cast<MeshLib::Mesh&>(mesh), "NodalForcesJump" + std::to_string(ig+1),
                    MeshLib::MeshItemType::Node, GlobalDim);
            assert(_process_data.vec_mesh_prop_nodal_forces_jump[ig]->size() ==
                GlobalDim * mesh.getNumberOfNodes());
        }

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
}

template <int GlobalDim>
void ThermoHydroMechanicsProcess<GlobalDim>::computeSecondaryVariableConcrete(
    const double t, GlobalVector const& x, int const process_id)
{
    DBUG("Compute the secondary variables for ThermoHydroMechanicsProcess.");
    const auto& dof_table = getDOFTable(process_id);

    {
        ProcessLib::ProcessVariable const& pv =
            getProcessVariables(process_id)[0];

        GlobalExecutor::executeSelectedMemberOnDereferenced(
            &ThermoHydroMechanicsLocalAssemblerInterface::computeSecondaryVariable,
            _local_assemblers, pv.getActiveElementIDs(),
            dof_table, t, x, _coupled_solutions);
    }

    // Copy displacement jumps in a solution vector to mesh property
    // Remark: the copy is required because mesh properties for primary
    // variables are set during output and are not ready yet when this function
    // is called.
    auto const n_enrich_vars = _vec_fracture_nodes.size() + _vec_junction_nodes.size();
    std::vector<int> vec_g_variable_id(n_enrich_vars);
    {
        const int monolithic_process_id = 0;
        auto const& pvs = getProcessVariables(monolithic_process_id);
        for (unsigned i=0; i<n_enrich_vars; i++)
        {
            auto const displacement_jump_name_i = "displacement_jump" + std::to_string(i+1);
            auto const it =
                std::find_if(pvs.begin(), pvs.end(), [&](ProcessVariable const& pv) {
                    return pv.getName() == displacement_jump_name_i;
                });
            if (it == pvs.end())
            {
                OGS_FATAL(
                    "Didn't find expected 'displacement_jump1' process "
                    "variable.");
            }
            vec_g_variable_id[i] = static_cast<int>(std::distance(pvs.begin(), it));
        }
    }

    MathLib::LinAlg::setLocalAccessibleVector(x);

    const int monolithic_process_id = 0;
    for (auto const g_variable_id : vec_g_variable_id)
    {
        ProcessVariable& pv_g =
            this->getProcessVariables(monolithic_process_id)[g_variable_id];
        auto const num_comp = pv_g.getNumberOfComponents();
        auto& mesh_prop_g = *MeshLib::getOrCreateMeshProperty<double>(
            _mesh, pv_g.getName(), MeshLib::MeshItemType::Node, num_comp);
        for (int component_id = 0; component_id < num_comp; ++component_id)
        {
            auto const& mesh_subset = dof_table.getMeshSubset(
                g_variable_id, component_id);
            auto const mesh_id = mesh_subset.getMeshID();
            for (auto const* node : mesh_subset.getNodes())
            {
                MeshLib::Location const l(mesh_id, MeshLib::MeshItemType::Node,
                                        node->getID());

                auto const global_index =
                    dof_table.getGlobalIndex(l, g_variable_id, component_id);
                mesh_prop_g[node->getID() * num_comp + component_id] =
                    x[global_index];
            }
        }
    }

    // compute nodal w and aperture
    MeshLib::PropertyVector<double>& vec_w = *_process_data.mesh_prop_nodal_w;
    MeshLib::PropertyVector<double>& vec_b = *_process_data.mesh_prop_nodal_b;

    Eigen::VectorXd g(GlobalDim);
    Eigen::VectorXd w(GlobalDim);
    for (unsigned frac_id=0; frac_id<_vec_fracture_nodes.size(); frac_id++)
    {
        auto const& R = _process_data.fracture_properties[frac_id]->R;
        auto const& b0 = _process_data.fracture_properties[frac_id]->aperture0;
        auto compute_nodal_aperture = [&](std::size_t const node_id,
                                        double const w_n) {
            // skip aperture computation for element-wise defined b0 because there
            // are jumps on the nodes between the element's values.
            if (dynamic_cast<ParameterLib::MeshElementParameter<double> const*>(
                    &b0))
            {
                return std::numeric_limits<double>::quiet_NaN();
            }

            ParameterLib::SpatialPosition x;
            x.setNodeID(node_id);
            return w_n + b0(/*time independent*/ 0, x)[0];
        };

        auto const& frac_elements = _vec_fracture_elements[frac_id];
        auto const& frac_nodes = _vec_fracture_nodes[frac_id];
        for (MeshLib::Node const* node : frac_nodes)
        {
            auto const node_id = node->getID();
            //---------------------------------------------------------------
            // find connecting fracture elements
            //---------------------------------------------------------------
            MeshLib::Element const* e = nullptr;
            for (unsigned i=0; i<node->getNumberOfElements(); i++)
            {
                if (node->getElement(i)->getDimension()==GlobalDim-1)
                {
                    e = node->getElement(i);
                    if (frac_elements.end() != std::find(frac_elements.begin(), frac_elements.end(), e))
                        break;
                }
            }
            assert(e!=nullptr);

            //---------------------------------------------------------------
            // calc levelsets
            //---------------------------------------------------------------
            Eigen::Vector3d const pt(e->getCenterOfGravity().getCoords());
            std::vector<FractureProperty*> e_fracture_props;
            std::unordered_map<int, int> e_fracID_to_local;
            unsigned tmpi = 0;
            for (auto fid :
                _process_data.vec_ele_connected_fractureIDs[e->getID()])
            {
                e_fracture_props.push_back(&*_process_data.fracture_properties[fid]);
                e_fracID_to_local.insert({fid, tmpi++});
            }
            std::vector<JunctionProperty*> e_junction_props;
            std::unordered_map<int, int> e_juncID_to_local;
            tmpi = 0;
            for (auto fid :
                _process_data.vec_ele_connected_junctionIDs[e->getID()])
            {
                e_junction_props.push_back(&_process_data.junction_properties[fid]);
                e_juncID_to_local.insert({fid, tmpi++});
            }
            std::vector<double> const levelsets(duGlobalEnrichments(
                frac_id, e_fracture_props, e_junction_props,
                e_fracID_to_local, pt));

            std::vector<int> localVarId_to_pvId(levelsets.size());
            for (unsigned i=0; i<e_fracture_props.size(); i++)
                localVarId_to_pvId[i] = e_fracture_props[i]->fracture_id;
            for (unsigned i=0; i<e_junction_props.size(); i++)
                localVarId_to_pvId[i+e_fracture_props.size()] = e_junction_props[i]->junction_id;

            //---------------------------------------------------------------
            // calc true g
            //---------------------------------------------------------------
            g.setZero();
            for (unsigned i = 0; i < levelsets.size(); i++)
            {
                auto const g_variable_id = vec_g_variable_id[localVarId_to_pvId[i]];
                ProcessVariable& pv_g =
                    this->getProcessVariables(monolithic_process_id)[g_variable_id];
                auto& mesh_prop_g = *MeshLib::getOrCreateMeshProperty<double>(
                    _mesh, pv_g.getName(), MeshLib::MeshItemType::Node, pv_g.getNumberOfComponents());
                for (int k = 0; k < GlobalDim; k++)
                {
                    g[k] += levelsets[i] * mesh_prop_g[node_id * GlobalDim + k];
                }
            }

            //---------------------------------------------------------------
            // calc new uperture
            //---------------------------------------------------------------
            w.noalias() = R * g;
            for (int k = 0; k < GlobalDim; k++)
            {
                vec_w[node_id * GlobalDim + k] = w[k];
            }

            vec_b[node_id] = compute_nodal_aperture(node_id, w[GlobalDim - 1]);
        }
    }
}

template <int GlobalDim>
bool ThermoHydroMechanicsProcess<GlobalDim>::isLinear() const
{
    return false;
}

template <int GlobalDim>
void ThermoHydroMechanicsProcess<GlobalDim>::assembleConcreteProcess(
    const double t, GlobalVector const& x, GlobalMatrix& M, GlobalMatrix& K,
    GlobalVector& b)
{
    DBUG("Assemble ThermoHydroMechanicsProcess.");

    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
        dof_table = {std::ref(*_local_to_global_index_map)};
    // Call global assembler for each local assembly item.
    GlobalExecutor::executeMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assemble, _local_assemblers,
        dof_table, t, x, M, K, b, _coupled_solutions);
}

template <int GlobalDim>
void ThermoHydroMechanicsProcess<GlobalDim>::assembleWithJacobianConcreteProcess(
    const double t, GlobalVector const& x, GlobalVector const& xdot,
    const double dxdot_dx, const double dx_dx, GlobalMatrix& M, GlobalMatrix& K,
    GlobalVector& b, GlobalMatrix& Jac)
{
    DBUG("AssembleWithJacobian ThermoHydroMechanicsProcess.");

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
    copyRhs(2, *_process_data.mesh_prop_nodal_forces);
    const auto n_pv_g = _process_data.fracture_properties.size() + _process_data.junction_properties.size();
    for (unsigned i=0; i<n_pv_g; i++)
        copyRhs(i+3, *_process_data.vec_mesh_prop_nodal_forces_jump[i]);
}

template <int GlobalDim>
void ThermoHydroMechanicsProcess<GlobalDim>::preTimestepConcreteProcess(
    GlobalVector const& x, double const t, double const dt,
    const int process_id)
{
    DBUG("PreTimestep ThermoHydroMechanicsProcess.");

    _process_data.dt = dt;
    _process_data.t = t;

    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];

    GlobalExecutor::executeSelectedMemberOnDereferenced(
        &ThermoHydroMechanicsLocalAssemblerInterface::preTimestep, _local_assemblers,
        pv.getActiveElementIDs(), *_local_to_global_index_map,
        x, t, dt);
}

// ------------------------------------------------------------------------------------
// template instantiation
// ------------------------------------------------------------------------------------
template class ThermoHydroMechanicsProcess<2>;
template class ThermoHydroMechanicsProcess<3>;

}  // namespace ThermoHydroMechanics
}  // namespace LIE
}  // namespace ProcessLib
