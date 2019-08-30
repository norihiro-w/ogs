
#include "CreateTHProcess.h"

#include <cassert>

#include "MaterialLib/Fluid/FluidPropertyHeaders.h"
#include "MaterialLib/Fluid/FluidProperties/CreateFluidProperties.h"
#include "MaterialLib/Fluid/FluidProperties/FluidProperties.h"
#include "MaterialLib/FractureModels/CreateCohesiveZoneModeI.h"
#include "MaterialLib/FractureModels/CreateLinearElasticIsotropic.h"
#include "MaterialLib/FractureModels/CreateMohrCoulomb.h"
#include "MaterialLib/FractureModels/Permeability/CreatePermeabilityModel.h"
#include "MaterialLib/SolidModels/CreateConstitutiveRelation.h"
#include "ParameterLib/Utils.h"
#include "ProcessLib/Output/CreateSecondaryVariables.h"

#include "THProcess.h"
#include "THProcessData.h"

namespace ProcessLib
{
namespace LIE
{
namespace TH
{
template <unsigned GlobalDim>
std::unique_ptr<Process> createTHProcess(
    std::string name,
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    boost::optional<ParameterLib::CoordinateSystem> const&
        /*local_coordinate_system*/,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{prj__processes__process__type}
    config.checkConfigParameter("type", "TH_WITH_LIE");
    DBUG("Create THProcess.");
    auto const staggered_scheme =
        //! \ogs_file_param{prj__processes__process__TH__coupling_scheme}
        config.getConfigParameterOptional<std::string>("coupling_scheme");
    const bool use_monolithic_scheme =
        !(staggered_scheme && (*staggered_scheme == "staggered"));

    // Process variables
    //! \ogs_file_param{prj__processes__process__TH__process_variables}
    auto const pv_conf = config.getConfigSubtree("process_variables");
    auto range =
        //! \ogs_file_param{prj__processes__process__TH__process_variables__process_variable}
        pv_conf.getConfigParameterList<std::string>("process_variable");
    std::vector<std::reference_wrapper<ProcessVariable>> p_T_process_variables;
    std::vector<std::reference_wrapper<ProcessVariable>> p_process_variables;
    std::vector<std::reference_wrapper<ProcessVariable>> T_process_variables;
    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>
        process_variables;

    unsigned n_var_du = 0;
    for (std::string const& pv_name : range)
    {
        if (pv_name != "pressure" && pv_name != "temperature")
        {
            OGS_FATAL(
                "Found a process variable name '%s'. It should be "
                "'displacement' or 'displacement_jumpN' or 'pressure' or 'temperature'");
        }

        if (pv_name.find("displacement_jump") == 0)
            n_var_du++;

        auto variable = std::find_if(variables.cbegin(), variables.cend(),
                                     [&pv_name](ProcessVariable const& v) {
                                         return v.getName() == pv_name;
                                     });

        if (variable == variables.end())
        {
            OGS_FATAL(
                "Could not find process variable '%s' in the provided "
                "variables "
                "list for config tag <%s>.",
                pv_name.c_str(), "process_variable");
        }
        DBUG("Found process variable '%s' for config tag <%s>.",
             variable->getName().c_str(), "process_variable");

        if (pv_name.find("displacement") != std::string::npos &&
            variable->getNumberOfComponents() != GlobalDim)
        {
            OGS_FATAL(
                "Number of components of the process variable '%s' is "
                "different "
                "from the displacement dimension: got %d, expected %d",
                variable->getName().c_str(),
                variable->getNumberOfComponents(),
                GlobalDim);
        }

        if (!use_monolithic_scheme)
        {
            if (pv_name == "pressure")
            {
                p_process_variables.emplace_back(
                    const_cast<ProcessVariable&>(*variable));
            }
            else if (pv_name == "temperature")
            {
                T_process_variables.emplace_back(
                    const_cast<ProcessVariable&>(*variable));
            }
        }
        else
        {
            p_T_process_variables.emplace_back(
                const_cast<ProcessVariable&>(*variable));
        }
    }

    if (!use_monolithic_scheme)
    {
        process_variables.push_back(std::move(p_process_variables));
        process_variables.push_back(std::move(T_process_variables));
    }
    else
    {
        process_variables.push_back(std::move(p_T_process_variables));
    }

    // Fluid properties
    auto const& fluid_config = config.getConfigSubtree("fluid");
    auto fluid_props =
        MaterialLib::Fluid::createFluidProperties(fluid_config);

    // auto solid_constitutive_relations =
    //     MaterialLib::Solids::createConstitutiveRelations<GlobalDim>(
    //         parameters, local_coordinate_system, config);

    // Intrinsic permeability
    auto& intrinsic_permeability = ParameterLib::findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__TH__intrinsic_permeability}
        "intrinsic_permeability", parameters, 1, &mesh);

    DBUG("Use '%s' as intrinsic permeability parameter.",
         intrinsic_permeability.name.c_str());

    // Storage coefficient
    auto& specific_storage = ParameterLib::findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__TH__specific_storage}
        "specific_storage", parameters, 1, &mesh);

    DBUG("Use '%s' as specific storage parameter.",
         specific_storage.name.c_str());

    // Biot coefficient
    auto& biot_coefficient = ParameterLib::findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__TH__biot_coefficient}
        "biot_coefficient", parameters, 1, &mesh);
    DBUG("Use '%s' as Biot coefficient parameter.",
         biot_coefficient.name.c_str());

    // Porosity
    auto& porosity = ParameterLib::findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__TH__porosity}
        "porosity", parameters, 1, &mesh);
    DBUG("Use '%s' as porosity parameter.", porosity.name.c_str());

    // Solid density
    auto& solid_density = ParameterLib::findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__TH__solid_density}
        "solid_density", parameters, 1, &mesh);
    DBUG("Use '%s' as solid density parameter.", solid_density.name.c_str());

    // linear thermal expansion coefficient for solid
    auto const& solid_linear_thermal_expansion_coefficient =
        ParameterLib::findParameter<double>(
            config,
            //! \ogs_file_param_special{prj__processes__process__THERMO_HYDRO_MECHANICS__solid_linear_thermal_expansion_coefficient}
            "solid_linear_thermal_expansion_coefficient", parameters, 1, &mesh);
    DBUG("Use '%s' as solid linear thermal expansion coefficient parameter.",
         solid_linear_thermal_expansion_coefficient.name.c_str());

    // specific heat capacity for solid
    auto& solid_specific_heat_capacity = ParameterLib::findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__THERMO_HYDRO_MECHANICS__solid_specific_heat_capacity}
        "solid_specific_heat_capacity", parameters, 1, &mesh);
    DBUG("Use '%s' as solid specific heat capacity parameter.",
         solid_specific_heat_capacity.name.c_str());

    // thermal conductivity for solid // currently only considers isotropic
    auto& solid_thermal_conductivity = ParameterLib::findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__THERMO_HYDRO_MECHANICS__solid_thermal_conductivity}
        "solid_thermal_conductivity", parameters, 1, &mesh);
    DBUG("Use '%s' as solid thermal conductivity parameter.",
         solid_thermal_conductivity.name.c_str());

    // Specific body force
    Eigen::Matrix<double, GlobalDim, 1> specific_body_force;
    {
        std::vector<double> const b =
            //! \ogs_file_param{prj__processes__process__TH__specific_body_force}
            config.getConfigParameter<std::vector<double>>(
                "specific_body_force");
        if (b.size() != GlobalDim)
        {
            OGS_FATAL(
                "The size of the specific body force vector does not match the "
                "displacement dimension. Vector size is %d, displacement "
                "dimension is %d",
                b.size(), GlobalDim);
        }

        std::copy_n(b.data(), b.size(), specific_body_force.data());
    }

    // Fracture constitutive relation.
    std::unique_ptr<MaterialLib::Fracture::FractureModelBase<GlobalDim>>
        fracture_model = nullptr;
    auto const opt_fracture_model_config =
        //! \ogs_file_param{prj__processes__process__TH__fracture_model}
        config.getConfigSubtreeOptional("fracture_model");
    // if (opt_fracture_model_config)
    // {
    //     auto& fracture_model_config = *opt_fracture_model_config;

    //     auto const frac_type =
    //         //! \ogs_file_param{prj__processes__process__TH__fracture_model__type}
    //         fracture_model_config.peekConfigParameter<std::string>("type");

    //     if (frac_type == "LinearElasticIsotropic")
    //     {
    //         fracture_model =
    //             MaterialLib::Fracture::createLinearElasticIsotropic<GlobalDim>(
    //                 parameters, fracture_model_config);
    //     }
    //     else if (frac_type == "MohrCoulomb")
    //     {
    //         fracture_model =
    //             MaterialLib::Fracture::createMohrCoulomb<GlobalDim>(
    //                 parameters, fracture_model_config);
    //     }
    //     else if (frac_type == "CohesiveZoneModeI")
    //     {
    //         fracture_model = MaterialLib::Fracture::CohesiveZoneModeI::
    //             createCohesiveZoneModeI<GlobalDim>(parameters,
    //                                                fracture_model_config);
    //     }
    //     else
    //     {
    //         OGS_FATAL(
    //             "Cannot construct fracture constitutive relation of given type "
    //             "'%s'.",
    //             frac_type.c_str());
    //     }
    // }

    // Fracture properties
    std::vector<std::unique_ptr<FractureProperty>> fracture_properties;
    for (auto fracture_properties_config :
        //! \ogs_file_param{prj__processes__process__TH__fracture_properties}
        config.getConfigSubtreeList("fracture_properties"))
    {
        auto frac_prop = std::make_unique<ProcessLib::LIE::FractureProperty>(
            fracture_properties.size(),
            //! \ogs_file_param{prj__processes__process__TH__fracture_properties__material_id}
            fracture_properties_config.getConfigParameter<int>("material_id"),
            ParameterLib::findParameter<double>(
                //! \ogs_file_param_special{prj__processes__process__TH__fracture_properties__initial_aperture}
                fracture_properties_config, "initial_aperture", parameters, 1,
                &mesh),
            ParameterLib::findParameter<double>(
                //! \ogs_file_param_special{prj__processes__process__TH__fracture_properties__specific_storage}
                fracture_properties_config, "specific_storage", parameters, 1,
                &mesh),
            ParameterLib::findParameter<double>(
                //! \ogs_file_param_special{prj__processes__process__TH__fracture_properties__biot_coefficient}
                fracture_properties_config, "biot_coefficient", parameters, 1,
                &mesh),
            ParameterLib::findParameter<double>(
                //! \ogs_file_param_special{prj__processes__process__TH__fracture_properties__initial_fracture_effective_stress}
                fracture_properties_config, "initial_fracture_effective_stress", parameters, GlobalDim,
                &mesh));
        if (frac_prop->aperture0.isTimeDependent())
        {
            OGS_FATAL(
                "The initial aperture parameter '%s' must not be "
                "time-dependent.",
                frac_prop->aperture0.name.c_str());
        }

        auto permeability_model_config =
            //! \ogs_file_param{prj__processes__process__TH__fracture_properties__permeability_model}
            fracture_properties_config.getConfigSubtree("permeability_model");
        frac_prop->permeability_model =
            MaterialLib::Fracture::Permeability::createPermeabilityModel(
                permeability_model_config);

        fracture_properties.emplace_back(std::move(frac_prop));
    }

    // initial effective stress in matrix
    auto& initial_effective_stress = ParameterLib::findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__TH__initial_effective_stress}
        "initial_effective_stress", parameters,
        GlobalDim==2 ? 4 : 6, &mesh);
    DBUG("Use '%s' as initial effective stress parameter.",
         initial_effective_stress.name.c_str());

    // deactivation of matrix elements in flow
    auto opt_deactivate_matrix_in_flow =
        //! \ogs_file_param{prj__processes__process__TH__deactivate_matrix_in_flow}
        config.getConfigParameterOptional<bool>("deactivate_matrix_in_flow");
    bool const deactivate_matrix_in_flow =
        opt_deactivate_matrix_in_flow && *opt_deactivate_matrix_in_flow;
    ;
    if (deactivate_matrix_in_flow)
        INFO("Deactivate matrix elements in flow calculation.");

    THProcessData<GlobalDim> process_data{
        materialIDs(mesh),
        std::move(fluid_props),
        solid_density,
        solid_linear_thermal_expansion_coefficient,
        solid_specific_heat_capacity,
        solid_thermal_conductivity,
        intrinsic_permeability,
        specific_storage,
        biot_coefficient,
        porosity,
        specific_body_force,
        std::move(fracture_model),
        std::move(fracture_properties),
        initial_effective_stress,
        deactivate_matrix_in_flow};

    SecondaryVariableCollection secondary_variables;

    NumLib::NamedFunctionCaller named_function_caller(
        {"TH"});

    ProcessLib::createSecondaryVariables(config, secondary_variables,
                                         named_function_caller);

    return std::make_unique<THProcess<GlobalDim>>(
        std::move(name), mesh, std::move(jacobian_assembler), parameters,
        integration_order, std::move(process_variables),
        std::move(process_data), std::move(secondary_variables),
        std::move(named_function_caller), use_monolithic_scheme);
}

template std::unique_ptr<Process> createTHProcess<1>(
    std::string name,
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    boost::optional<ParameterLib::CoordinateSystem> const&
        local_coordinate_system,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config);
template std::unique_ptr<Process> createTHProcess<2>(
    std::string name,
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    boost::optional<ParameterLib::CoordinateSystem> const&
        local_coordinate_system,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config);
template std::unique_ptr<Process> createTHProcess<3>(
    std::string name,
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    boost::optional<ParameterLib::CoordinateSystem> const&
        local_coordinate_system,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config);

}  // namespace TH
}  // namespace LIE
}  // namespace ProcessLib
