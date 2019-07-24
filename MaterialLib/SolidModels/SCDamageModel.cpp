/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "SCDamageModel.h"

#include <iostream>
#include <limits>

#include "BaseLib/Error.h"
#include "MaterialLib/PhysicalConstant.h"

namespace MaterialLib
{
namespace Solids
{
namespace SCDamage
{

namespace detail
{
/// Lamé's first parameter.
double lambda(double E, double nu)
{
    return E * nu / (1 + nu) / (1 - 2 * nu);
}

/// Lamé's second parameter, the shear modulus.
double mu(double E, double nu)
{
    return E / (2 * (1 + nu));
}

double bulk_modulus(double E, double nu)
{
    return E / (3 * (1 - 2 * nu));
}

}

template <int DisplacementDim>
boost::optional<std::tuple<typename SCDamageModel<DisplacementDim>::KelvinVector,
                           std::unique_ptr<typename MechanicsBase<
                               DisplacementDim>::MaterialStateVariables>,
                           typename SCDamageModel<DisplacementDim>::KelvinMatrix>>
SCDamageModel<DisplacementDim>::integrateStress(
    double const t, ProcessLib::SpatialPosition const& x, double const /*dt*/,
    KelvinVector const& eps_prev, KelvinVector const& eps,
    KelvinVector const& sigma_prev,
    typename MechanicsBase<DisplacementDim>::MaterialStateVariables const&
    /*material_state_variables*/,
    double const T) const
{
    //using Invariants = MathLib::KelvinVector::Invariants<KelvinVectorSize>;

    const auto damaged = _damage_state(t, x)[0];
    const auto C = this->getElasticTensor(t, x, T, damaged);
    KelvinVector sigma_try = sigma_prev + C * (eps - eps_prev);

    return {std::make_tuple(sigma_try, createMaterialStateVariables(),
                            C)};
}


template <int DisplacementDim>
typename SCDamageModel<DisplacementDim>::KelvinMatrix
SCDamageModel<DisplacementDim>::getElasticTensor(
    double const t, ProcessLib::SpatialPosition const& x,
    double const T, double const damaged) const
{
    ProcessLib::ParameterArguments args;
    args.t = t;
    args.pos = x;
    args.exargs.emplace("T", T);
    auto const E = damaged==0 ? _mp._youngs_modulus(args)[0] : _mp._youngs_modulus_damaged(args)[0];
    auto const nu = _mp._poissons_ratio(args)[0];
    return elasticTangentStiffness<DisplacementDim>(detail::lambda(E, nu),
                                                    detail::mu(E, nu));
}

template class SCDamageModel<2>;
template class SCDamageModel<3>;

}
}  // namespace Solids
}  // namespace MaterialLib
