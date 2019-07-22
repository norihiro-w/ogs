/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "ProcessLib/Utils/ProcessUtils.h"  // required for findParameter
#include "SCDamageModel.h"

namespace MaterialLib
{
namespace Solids
{
namespace SCDamage
{

template <int DisplacementDim>
std::unique_ptr<SCDamageModel<DisplacementDim>>
createSCDamageModel(
    std::vector<std::unique_ptr<ProcessLib::ParameterBase>> const& parameters,
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{material__solid__constitutive_relation__type}
    config.checkConfigParameter("type", "SCDamageModel");
    DBUG("Create SCDamageModel material");

    // Youngs modulus
    auto& youngs_modulus = ProcessLib::findParameter<double>(
        //! \ogs_file_param_special{material__solid__constitutive_relation__SCDamageModel__youngs_modulus}
        config, "youngs_modulus", parameters, 1);

    DBUG("Use '%s' as youngs_modulus parameter.", youngs_modulus.name.c_str());

    // Poissons ratio
    auto& poissons_ratio = ProcessLib::findParameter<double>(
        //! \ogs_file_param_special{material__solid__constitutive_relation__SCDamageModel__poissons_ratio}
        config, "poissons_ratio", parameters, 1);

    DBUG("Use '%s' as poissons_ratio parameter.", poissons_ratio.name.c_str());

    // Youngs modulus of damaged
    auto& youngs_modulus_damaged = ProcessLib::findParameter<double>(
        //! \ogs_file_param_special{material__solid__constitutive_relation__SCDamageModel__youngs_modulus_damaged}
        config, "youngs_modulus_damaged", parameters, 1);

    DBUG("Use '%s' as youngs_modulus_damaged parameter.", youngs_modulus_damaged.name.c_str());

    // Fricition angle
    auto& friction_angle = ProcessLib::findParameter<double>(
        //! \ogs_file_param_special{material__solid__constitutive_relation__SCDamageModel__friction_angle}
        config, "friction_angle", parameters, 1);

    // Cohesion
    auto& cohesion = ProcessLib::findParameter<double>(
        //! \ogs_file_param_special{material__solid__constitutive_relation__SCDamageModel__cohesion}
        config, "cohesion", parameters, 1);

    DBUG("Use '%s' as youngs_modulus_damaged parameter.", youngs_modulus_damaged.name.c_str());

    // Damage state
    auto& damage_state = ProcessLib::findParameter<double>(
        //! \ogs_file_param_special{material__solid__constitutive_relation__SCDamageModel__damage_state}
        config, "damage_state", parameters, 1);

    DBUG("Use '%s' as damage_state parameter.", damage_state.name.c_str());

    typename SCDamageModel<DisplacementDim>::MaterialProperties mp{
        youngs_modulus, youngs_modulus_damaged, poissons_ratio,
        friction_angle, cohesion};

    return std::make_unique<SCDamageModel<DisplacementDim>>(mp, damage_state);
}

}
}  // namespace Solids
}  // namespace MaterialLib
