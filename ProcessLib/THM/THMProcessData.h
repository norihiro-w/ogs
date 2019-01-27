/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <memory>
#include <utility>

#include <Eigen/Dense>

#include "MaterialLib/Fluid/FluidProperties/FluidProperties.h"
#include "MaterialLib/Fluid/FluidPropertyHeaders.h"
#include "ProcessLib/Parameter/Parameter.h"

namespace MaterialLib
{
namespace Solids
{
template <int DisplacementDim>
struct MechanicsBase;
}
}  // namespace MaterialLib
namespace ProcessLib
{
namespace THM
{
template <int DisplacementDim>
struct THMProcessData
{
    THMProcessData(
        MeshLib::PropertyVector<int> const* const material_ids_,
        std::map<int,
                 std::unique_ptr<
                     MaterialLib::Solids::MechanicsBase<DisplacementDim>>>&&
            solid_materials_,
        std::unique_ptr<MaterialLib::Fluid::FluidProperties>&& fluid_props_,
        Parameter<double> const& solid_density_,
        Parameter<double> const& solid_linear_thermal_expansion_coefficient_,
        Parameter<double> const& solid_specific_heat_capacity_,
        Parameter<double> const& solid_thermal_conductivity_,
        Parameter<double> const& intrinsic_permeability_,
        Parameter<double> const& specific_storage_,
        Parameter<double> const& biot_coefficient_,
        Parameter<double> const& porosity_,
        Eigen::Matrix<double, DisplacementDim, 1>
            specific_body_force_,
        bool const reset_strain_)
        // double const reference_temperature_)
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
          reset_strain(reset_strain_)
    {
    }

    THMProcessData(THMProcessData&& other) =
        default;

    //! Copies are forbidden.
    THMProcessData(THMProcessData const&) =
        delete;

    //! Assignments are not needed.
    void operator=(THMProcessData const&) = delete;

    //! Assignments are not needed.
    void operator=(THMProcessData&&) = delete;

    MeshLib::PropertyVector<int> const* const material_ids;

    /// The constitutive relation for the mechanical part.
    /// \note Linear elasticity is the only supported one in the moment.
    std::map<
        int,
        std::unique_ptr<MaterialLib::Solids::MechanicsBase<DisplacementDim>>>
        solid_materials;

    /// Fluid properties
    std::unique_ptr<MaterialLib::Fluid::FluidProperties> fluid_props;
    MaterialLib::Fluid::FluidProperty const& fluid_density;
    MaterialLib::Fluid::FluidProperty const& fluid_viscosity;
    MaterialLib::Fluid::FluidProperty const& fluid_specific_heat_capacity;
    MaterialLib::Fluid::FluidProperty const& fluid_thermal_conductivity;

    /// Solid properties
    Parameter<double> const& solid_density;
    Parameter<double> const& solid_linear_thermal_expansion_coefficient;
    Parameter<double> const& solid_specific_heat_capacity;
    Parameter<double> const& solid_thermal_conductivity;

    /// Permeability of the solid. A scalar quantity, Parameter<double>.
    Parameter<double> const& intrinsic_permeability;
    /// Volumetric average specific storage of the solid and fluid phases.
    /// A scalar quantity, Parameter<double>.
    Parameter<double> const& specific_storage;
    /// Biot coefficient. A scalar quantity, Parameter<double>.
    Parameter<double> const& biot_coefficient;
    /// Porosity of the solid. A scalar quantity, Parameter<double>.
    Parameter<double> const& porosity;

    /// Specific body forces applied to solid and fluid.
    /// It is usually used to apply gravitational forces.
    /// A vector of displacement dimension's length.
    Eigen::Matrix<double, DisplacementDim, 1> const specific_body_force;

    bool const reset_strain;

    double dt = 0.0;
    double t = 0.0;

    // double const reference_temperature;

    // MeshLib::PropertyVector<double>* pressure_interpolated = nullptr;
    // MeshLib::PropertyVector<double>* temperature_interpolated = nullptr;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

}  // namespace THM
}  // namespace ProcessLib
