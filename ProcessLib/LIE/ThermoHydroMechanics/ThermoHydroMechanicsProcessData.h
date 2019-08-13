/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <memory>
#include <utility>

#include <Eigen/Dense>

#include "MeshLib/ElementStatus.h"
#include "MeshLib/PropertyVector.h"

#include "MaterialLib/Fluid/FluidProperties/FluidProperties.h"
#include "MaterialLib/Fluid/FluidPropertyHeaders.h"
#include "MaterialLib/FractureModels/FractureModelBase.h"
#include "MaterialLib/SolidModels/MechanicsBase.h"

#include "ProcessLib/LIE/Common/FractureProperty.h"
#include "ProcessLib/LIE/Common/JunctionProperty.h"

namespace MeshLib
{
class Element;
}

namespace ProcessLib
{
namespace LIE
{
namespace ThermoHydroMechanics
{
template <unsigned GlobalDim>
struct ThermoHydroMechanicsProcessData
{
    ThermoHydroMechanicsProcessData(
        MeshLib::PropertyVector<int> const* const material_ids_,
        std::map<
            int,
            std::unique_ptr<MaterialLib::Solids::MechanicsBase<GlobalDim>>>&&
            solid_materials_,
        std::unique_ptr<MaterialLib::Fluid::FluidProperties>&& fluid_props_,
        ParameterLib::Parameter<double> const& solid_density_,
        ParameterLib::Parameter<double> const& solid_linear_thermal_expansion_coefficient_,
        ParameterLib::Parameter<double> const& solid_specific_heat_capacity_,
        ParameterLib::Parameter<double> const& solid_thermal_conductivity_,
        ParameterLib::Parameter<double> const& intrinsic_permeability_,
        ParameterLib::Parameter<double> const& specific_storage_,
        ParameterLib::Parameter<double> const& biot_coefficient_,
        ParameterLib::Parameter<double> const& porosity_,
        Eigen::Matrix<double, GlobalDim, 1>
            specific_body_force_,
        std::unique_ptr<MaterialLib::Fracture::FractureModelBase<GlobalDim>>&&
            fracture_model,
        std::vector<std::unique_ptr<FractureProperty>>&& fracture_properties,
        ParameterLib::Parameter<double> const& initial_effective_stress_,
        bool const deactivate_matrix_in_flow_)
        : material_ids(material_ids_),
          solid_materials{std::move(solid_materials_)},
          fluid_props{std::move(fluid_props_)},
          fluid_density(fluid_props->getProperty(MaterialLib::Fluid::FluidPropertyType::Density)),
          fluid_viscosity(fluid_props->getProperty(MaterialLib::Fluid::FluidPropertyType::Viscosity)),
          fluid_specific_heat_capacity(fluid_props->getProperty(MaterialLib::Fluid::FluidPropertyType::HeatCapacity)),
          fluid_thermal_conductivity(fluid_props->getProperty(MaterialLib::Fluid::FluidPropertyType::ThermalConductivity)),
          solid_density(solid_density_),
          solid_linear_thermal_expansion_coefficient(
              solid_linear_thermal_expansion_coefficient_),
          solid_specific_heat_capacity(solid_specific_heat_capacity_),
          solid_thermal_conductivity(solid_thermal_conductivity_),
          intrinsic_permeability(intrinsic_permeability_),
          specific_storage(specific_storage_),
          biot_coefficient(biot_coefficient_),
          porosity(porosity_),
          specific_body_force(std::move(specific_body_force_)),
          fracture_model{std::move(fracture_model)},
          fracture_properties{std::move(fracture_properties)},
          initial_effective_stress(initial_effective_stress_),
          deactivate_matrix_in_flow(deactivate_matrix_in_flow_)
    {
    }

    ThermoHydroMechanicsProcessData(ThermoHydroMechanicsProcessData&& other) = default;

    //! Copies are forbidden.
    ThermoHydroMechanicsProcessData(ThermoHydroMechanicsProcessData const&) = delete;

    //! Assignments are not needed.
    void operator=(ThermoHydroMechanicsProcessData const&) = delete;

    //! Assignments are not needed.
    void operator=(ThermoHydroMechanicsProcessData&&) = delete;

    /// Fluid properties
    std::unique_ptr<MaterialLib::Fluid::FluidProperties> fluid_props;
    MaterialLib::Fluid::FluidProperty const& fluid_density;
    MaterialLib::Fluid::FluidProperty const& fluid_viscosity;
    MaterialLib::Fluid::FluidProperty const& fluid_specific_heat_capacity;
    MaterialLib::Fluid::FluidProperty const& fluid_thermal_conductivity;

    /// Solid properties
    ParameterLib::Parameter<double> const& solid_density;
    ParameterLib::Parameter<double> const& solid_linear_thermal_expansion_coefficient;
    ParameterLib::Parameter<double> const& solid_specific_heat_capacity;
    ParameterLib::Parameter<double> const& solid_thermal_conductivity;

    /// Bulk properties
    MeshLib::PropertyVector<int> const* const material_ids;
    std::map<int,
             std::unique_ptr<MaterialLib::Solids::MechanicsBase<GlobalDim>>>
        solid_materials;
    ParameterLib::Parameter<double> const& intrinsic_permeability;
    ParameterLib::Parameter<double> const& specific_storage;
    ParameterLib::Parameter<double> const& biot_coefficient;
    ParameterLib::Parameter<double> const& porosity;
    Eigen::Matrix<double, GlobalDim, 1> const specific_body_force;

    /// Fracture properties
    std::unique_ptr<MaterialLib::Fracture::FractureModelBase<GlobalDim>>
        fracture_model;
    std::vector<std::unique_ptr<FractureProperty>> fracture_properties;
    std::vector<JunctionProperty> junction_properties;
    MeshLib::PropertyVector<int> const* mesh_prop_materialIDs = nullptr;
    std::vector<int> map_materialID_to_fractureID;
    // a table of connected fracture IDs for each element
    std::vector<std::vector<int>> vec_ele_connected_fractureIDs;
    std::vector<std::vector<int>> vec_ele_connected_junctionIDs;

    ParameterLib::Parameter<double> const& initial_effective_stress;

    bool const deactivate_matrix_in_flow;
    std::unique_ptr<MeshLib::ElementStatus> p_element_status;
    ParameterLib::Parameter<double> const* p0 = nullptr;

    double dt = 0.0;
    double t = 0.0;

    // mesh properties for output
    MeshLib::PropertyVector<double>* mesh_prop_stress_xx = nullptr;
    MeshLib::PropertyVector<double>* mesh_prop_stress_yy = nullptr;
    MeshLib::PropertyVector<double>* mesh_prop_stress_zz = nullptr;
    MeshLib::PropertyVector<double>* mesh_prop_stress_xy = nullptr;
    MeshLib::PropertyVector<double>* mesh_prop_stress_yz = nullptr;
    MeshLib::PropertyVector<double>* mesh_prop_stress_xz = nullptr;
    MeshLib::PropertyVector<double>* mesh_prop_strain_xx = nullptr;
    MeshLib::PropertyVector<double>* mesh_prop_strain_yy = nullptr;
    MeshLib::PropertyVector<double>* mesh_prop_strain_zz = nullptr;
    MeshLib::PropertyVector<double>* mesh_prop_strain_xy = nullptr;
    MeshLib::PropertyVector<double>* mesh_prop_strain_yz = nullptr;
    MeshLib::PropertyVector<double>* mesh_prop_strain_xz = nullptr;
    MeshLib::PropertyVector<double>* mesh_prop_velocity = nullptr;
    MeshLib::PropertyVector<double>* mesh_prop_b = nullptr;
    MeshLib::PropertyVector<double>* mesh_prop_k_f = nullptr;
    MeshLib::PropertyVector<double>* mesh_prop_w_n = nullptr;
    MeshLib::PropertyVector<double>* mesh_prop_w_s = nullptr;
    MeshLib::PropertyVector<double>* mesh_prop_w_s2 = nullptr;
    MeshLib::PropertyVector<double>* mesh_prop_fracture_stress_shear = nullptr;
    MeshLib::PropertyVector<double>* mesh_prop_fracture_stress_shear2 = nullptr;
    MeshLib::PropertyVector<double>* mesh_prop_fracture_stress_normal = nullptr;
    MeshLib::PropertyVector<double>* mesh_prop_fracture_shear_failure = nullptr;
    MeshLib::PropertyVector<double>* mesh_prop_nodal_w = nullptr;
    MeshLib::PropertyVector<double>* mesh_prop_nodal_b = nullptr;
    MeshLib::PropertyVector<double>* mesh_prop_nodal_p = nullptr;
    MeshLib::PropertyVector<double>* mesh_prop_nodal_T = nullptr;

    MeshLib::PropertyVector<double>* mesh_prop_nodal_forces = nullptr;
    std::vector<MeshLib::PropertyVector<double>*> vec_mesh_prop_nodal_forces_jump;
    MeshLib::PropertyVector<double>* mesh_prop_hydraulic_flow = nullptr;
    MeshLib::PropertyVector<double>* mesh_prop_thermal_flow = nullptr;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

}  // namespace ThermoHydroMechanics
}  // namespace LIE
}  // namespace ProcessLib
