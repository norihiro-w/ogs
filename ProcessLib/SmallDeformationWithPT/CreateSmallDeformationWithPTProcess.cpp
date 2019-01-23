/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateSmallDeformationWithPTProcess.h"

#include <cassert>

#include "MaterialLib/SolidModels/CreateConstitutiveRelation.h"
#include "MaterialLib/SolidModels/MechanicsBase.h"
#include "ProcessLib/Output/CreateSecondaryVariables.h"
#include "ProcessLib/Utils/ProcessUtils.h"

#include "SmallDeformationWithPTProcess.h"
#include "SmallDeformationWithPTProcessData.h"

namespace ProcessLib
{
namespace SmallDeformationWithPT
{
template <int DisplacementDim>
std::unique_ptr<Process> createSmallDeformationWithPTProcess(
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{prj__processes__process__type}
    config.checkConfigParameter("type", "SMALL_DEFORMATION_WITH_PT");
    DBUG("Create SmallDeformationWithPTProcess.");

    //  auto const staggered_scheme =
    //     //!
    //     \ogs_file_param{prj__processes__process__SMALL_DEFORMATION_WITH_PT__coupling_scheme}
    //     config.getConfigParameterOptional<std::string>("coupling_scheme");
    // const bool use_monolithic_scheme =
    //     !(staggered_scheme && (*staggered_scheme == "staggered"));

    // Process variable.

    //! \ogs_file_param{prj__processes__process__SMALL_DEFORMATION_WITH_PT__process_variables}
    auto const pv_config = config.getConfigSubtree("process_variables");

    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>
        process_variables;
    auto per_process_variables = findProcessVariables(
        variables, pv_config,
        {//! \ogs_file_param_special{prj__processes__process__SMALL_DEFORMATION_WITH_PT__process_variables__displacement}
         "displacement"});
    auto* variable_u = &per_process_variables[0].get();
    process_variables.push_back(std::move(per_process_variables));

    DBUG("Associate displacement with process variable '%s'.",
         variable_u->getName().c_str());

    if (variable_u->getNumberOfComponents() != DisplacementDim)
    {
        OGS_FATAL(
            "Number of components of the process variable '%s' is different "
            "from the displacement dimension: got %d, expected %d",
            variable_u->getName().c_str(),
            variable_u->getNumberOfComponents(),
            DisplacementDim);
    }

    //! \ogs_file_param{prj__processes__process__SMALL_DEFORMATION_WITH_PT__constitutive_relation}
    config.peekConfigParameter<std::string>("constitutive_relation");
    auto solid_constitutive_relations =
        MaterialLib::Solids::createConstitutiveRelations<DisplacementDim>(
            parameters, config);

    // pore pressure
    auto& p0 = findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__SMALL_DEFORMATION_WITH_PT__pressure_prev}
        "pressure_prev", parameters, 1);
    DBUG("Use '%s' as pressure_prev parameter.", p0.name.c_str());

    auto& p1 = findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__SMALL_DEFORMATION_WITH_PT__pressure}
        "pressure", parameters, 1);
    DBUG("Use '%s' as pressure parameter.", p1.name.c_str());

    // temperature
    auto& T0 = findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__SMALL_DEFORMATION_WITH_PT__temperature_prev}
        "temperature_prev", parameters, 1);
    DBUG("Use '%s' as temperature_prev parameter.", T0.name.c_str());

    auto& T1 = findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__SMALL_DEFORMATION_WITH_PT__temperature}
        "temperature", parameters, 1);
    DBUG("Use '%s' as temperature parameter.", T1.name.c_str());

    // solid density
    auto& solid_density = findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__SMALL_DEFORMATION_WITH_PT__solid_density}
        "solid_density", parameters, 1);
    DBUG("Use '%s' as solid density parameter.", solid_density.name.c_str());

    // Linear thermal expansion coefficient
    auto& solid_linear_thermal_expansion_coefficient = findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__SMALL_DEFORMATION_WITH_PT__solid_linear_thermal_expansion_coefficient}
        "solid_linear_thermal_expansion_coefficient", parameters, 1);
    DBUG("Use '%s' as linear thermal expansion coefficient.",
         solid_linear_thermal_expansion_coefficient.name.c_str());

    // Specific body force
    Eigen::Matrix<double, DisplacementDim, 1> specific_body_force;
    {
        std::vector<double> const b =
            //! \ogs_file_param{prj__processes__process__SMALL_DEFORMATION_WITH_PT__specific_body_force}
            config.getConfigParameter<std::vector<double>>(
                "specific_body_force");
        if (b.size() != DisplacementDim)
            OGS_FATAL(
                "The size of the specific body force vector does not match the "
                "displacement dimension. Vector size is %d, displacement "
                "dimension is %d",
                b.size(), DisplacementDim);

        std::copy_n(b.data(), b.size(), specific_body_force.data());
    }

    // biot coefficient
    auto& biot_coefficient = findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__SMALL_DEFORMATION_WITH_PT__biot_coefficient}
        "biot_coefficient", parameters, 1);
    DBUG("Use '%s' as biot coefficient parameter.",
         biot_coefficient.name.c_str());

    // biot coefficient
    auto& fluid_density = findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__SMALL_DEFORMATION_WITH_PT__fluid_density}
        "fluid_density", parameters, 1);
    DBUG("Use '%s' as fluid density parameter.", fluid_density.name.c_str());

    // porosity
    auto& porosity = findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__SMALL_DEFORMATION_WITH_PT__porosity}
        "porosity", parameters, 1);
    DBUG("Use '%s' as porosity parameter.", porosity.name.c_str());

    // external update
    //! \ogs_file_param{prj__processes__process__SMALL_DEFORMATION_WITH_PT__external_update}
    auto const config_import_mesh_properties =
        config.getConfigSubtreeOptional("import_mesh_properties");
    std::vector<std::pair<MeshLib::PropertyVector<double>*, std::string>>
        vec_import_properties;
    if (config_import_mesh_properties)
    {
        for (auto& config_property :
             config_import_mesh_properties->getConfigSubtreeList("property"))
        {
            auto const property_name =
                config_property.getConfigParameter<std::string>("name");
            DBUG("Importing property %s", property_name.c_str());
            auto const property_type =
                config_property.getConfigParameter<std::string>("type");

            MeshLib::PropertyVector<double>* property = nullptr;
            if (property_type == "MeshNode")
            {
                property = MeshLib::getOrCreateMeshProperty<double>(
                    const_cast<MeshLib::Mesh&>(mesh), property_name,
                    MeshLib::MeshItemType::Node, 1);
                property->resize(mesh.getNumberOfNodes());
            }
            else if (property_type == "MeshElement")
            {
                property = MeshLib::getOrCreateMeshProperty<double>(
                    const_cast<MeshLib::Mesh&>(mesh), property_name,
                    MeshLib::MeshItemType::Cell, 1);
                property->resize(mesh.getNumberOfElements());
            }

            //! \ogs_file_param{prj__processes__process__SMALL_DEFORMATION_WITH_PT__external_update__file}
            auto const file_name =
                config_property.getConfigParameter<std::string>("file");
            DBUG("Using file_name %s", file_name.c_str());
            auto const file_path =
                BaseLib::joinPaths(BaseLib::getProjectDirectory(), file_name);

            vec_import_properties.emplace_back(
                std::make_pair(property, file_path));
        }
    }

    // export
    //! \ogs_file_param{prj__processes__process__SMALL_DEFORMATION_WITH_PT__external_update}
    auto const config_export =
        config.getConfigSubtreeOptional("export_mesh_properties");
    std::vector<std::pair<MeshLib::PropertyVector<double>*, std::string>>
        vec_export_properties;
    if (config_export)
    {
        for (auto& config_property :
             config_export->getConfigSubtreeList("property"))
        {
            auto const property_name =
                config_property.getConfigParameter<std::string>("name");
            DBUG("Exporting property %s", property_name.c_str());

            auto const* property =
                mesh.getProperties().getPropertyVector<double>(property_name);
            if (!property)
                OGS_FATAL("Mesh property %s not found", property_name.c_str());

            auto const file_name =
                config_property.getConfigParameter<std::string>("file");
            DBUG("Using file_name %s", file_name.c_str());
            auto const file_path =
                BaseLib::joinPaths(BaseLib::getProjectDirectory(), file_name);

            vec_export_properties.emplace_back(std::make_pair(
                const_cast<MeshLib::PropertyVector<double>*>(property),
                file_path));
        }
    }

    SmallDeformationWithPTProcessData<DisplacementDim> process_data{
        materialIDs(mesh),
        T0,
        T1,
        p0,
        p1,
        std::move(solid_constitutive_relations),
        solid_density,
        solid_linear_thermal_expansion_coefficient,
        specific_body_force,
        biot_coefficient,
        fluid_density,
        porosity,
        vec_import_properties,
        vec_export_properties};

    SecondaryVariableCollection secondary_variables;

    NumLib::NamedFunctionCaller named_function_caller(
        {"SmallDeformationWithPT_displacement"});

    ProcessLib::createSecondaryVariables(config, secondary_variables,
                                         named_function_caller);

    return std::make_unique<SmallDeformationWithPTProcess<DisplacementDim>>(
        mesh, std::move(jacobian_assembler), parameters, integration_order,
        std::move(process_variables), std::move(process_data),
        std::move(secondary_variables), std::move(named_function_caller));
}

template std::unique_ptr<Process> createSmallDeformationWithPTProcess<2>(
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config);

template std::unique_ptr<Process> createSmallDeformationWithPTProcess<3>(
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config);

}  // namespace SmallDeformationWithPT
}  // namespace ProcessLib
