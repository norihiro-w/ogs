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

#include <Eigen/Eigen>

#include "ProcessLib/Parameter/Parameter.h"

namespace MaterialLib
{
namespace Solids
{
template <int DisplacementDim>
struct MechanicsBase;
}
}
namespace ProcessLib
{
namespace SmallDeformationWithPT
{
template <int DisplacementDim>
struct SmallDeformationWithPTProcessData
{
    SmallDeformationWithPTProcessData(
        MeshLib::PropertyVector<int> const* const material_ids_,
        Parameter<double> const& T0_,
        Parameter<double> const& T1_,
        Parameter<double> const& p0_,
        Parameter<double> const& p1_,
        std::map<int,
                 std::unique_ptr<
                     MaterialLib::Solids::MechanicsBase<DisplacementDim>>>&&
            solid_materials_,
        Parameter<double> const& solid_density_,
        Parameter<double> const& solid_linear_thermal_expansion_coefficient_,
        Eigen::Matrix<double, DisplacementDim, 1> const& specific_body_force_,
        Parameter<double> const& biot_coefficient_,
        Parameter<double> const& fluid_density_,
        Parameter<double> const& porosity_)
        : material_ids(material_ids_),
          T0(T0_), T1(T1_), p0(p0_), p1(p1_),
          solid_materials{std::move(solid_materials_)},
          solid_density(solid_density_),
          solid_linear_thermal_expansion_coefficient(
              solid_linear_thermal_expansion_coefficient_),
          specific_body_force(specific_body_force_),
          biot_coefficient(biot_coefficient_),
          fluid_density(fluid_density_),
          porosity(porosity_)
    {
    }

    SmallDeformationWithPTProcessData(SmallDeformationWithPTProcessData&& other) = default;

    //! Copies are forbidden.
    SmallDeformationWithPTProcessData(SmallDeformationWithPTProcessData const&) = delete;

    //! Assignments are not needed.
    void operator=(SmallDeformationWithPTProcessData const&) = delete;

    //! Assignments are not needed.
    void operator=(SmallDeformationWithPTProcessData&&) = delete;

    MeshLib::PropertyVector<int> const* const material_ids;

    Parameter<double> const& T0;
    Parameter<double> const& T1;
    Parameter<double> const& p0;
    Parameter<double> const& p1;

    /// The constitutive relation for the mechanical part.
    std::map<
        int,
        std::unique_ptr<MaterialLib::Solids::MechanicsBase<DisplacementDim>>>
        solid_materials;

    Parameter<double> const& solid_density;
    Parameter<double> const& solid_linear_thermal_expansion_coefficient;
    Eigen::Matrix<double, DisplacementDim, 1> const specific_body_force;
    Parameter<double> const& biot_coefficient;
    Parameter<double> const& fluid_density;
    Parameter<double> const& porosity;

    double dt = 0;
    double t = 0;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

}  // namespace SmallDeformationWithPT
}  // namespace ProcessLib
