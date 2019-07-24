/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "SCDamageModel2.h"

#include <iostream>
#include <limits>

#include "BaseLib/Error.h"
#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "MaterialLib/PhysicalConstant.h"

namespace MaterialLib
{
namespace Solids
{
namespace SCDamage2
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
struct PhysicalStressWithInvariants final
{
    static int const KelvinVectorSize =
        MathLib::KelvinVector::KelvinVectorDimensions<DisplacementDim>::value;
    using Invariants = MathLib::KelvinVector::Invariants<KelvinVectorSize>;
    using KelvinVector =
        MathLib::KelvinVector::KelvinVectorType<DisplacementDim>;

    explicit PhysicalStressWithInvariants(KelvinVector const& stress)
        : value{stress},
          D{Invariants::deviatoric_projection * stress},
          I_1{Invariants::trace(stress)},
          J_2{Invariants::J2(D)},
          J_3{Invariants::J3(D)}
    {
    }

    KelvinVector value;
    KelvinVector D;
    double I_1;
    double J_2;
    double J_3;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};


// template <int DisplacementDim>
// double yieldFunctionMC(MaterialProperties const& mp,
//                      PhysicalStressWithInvariants<DisplacementDim> const& s)
// {
//     assert(s.J_2 != 0);

//     double const phir = MathLib::to_radians(mp._friction_angle(t,x)[0]);
//     double const c = mp._cohesion(t,x)[0];

//     double theta = 1/3.*asin(-3.*std::sqrt(3.)*s.J_3/2./std::pow(s.J_2,1.5));
//     double m = cos(theta) - 1./std::sqrt(3.)*sin(theta)*sin(phir);
//     double beta = sin(phir)/m;
//     double kappa = cos(phir)/m*c;
//     double F = std::sqrt(J2) + beta * s.I_1/3. - kappa;

//     return F;
// }

template <int DisplacementDim>
double SCDamageModel2<DisplacementDim>::yieldFunctionMC(
    double const t, ProcessLib::SpatialPosition const& x, double const /*dt*/,
    KelvinVector const& sigma,
    typename MechanicsBase<DisplacementDim>::MaterialStateVariables const& /*material_state_variables*/
    ) const
{
//    using Invariants = MathLib::KelvinVector::Invariants<KelvinVectorSize>;
//    MaterialProperties const mp(t, x, _mp);
    auto& mp = _mp;

    PhysicalStressWithInvariants<DisplacementDim> const s{sigma};
    assert(s.J_2 != 0);

    double const phir = MathLib::to_radians(mp._friction_angle(t,x)[0]);
    double const c = mp._cohesion(t,x)[0];

    double sin_value = -3.*std::sqrt(3.)*s.J_3/2./std::pow(s.J_2,1.5);
    if (sin_value<-1. && -1-1e-6<sin_value)
        sin_value = -1;
    if (sin_value>1. && sin_value<1+1e-6)
        sin_value = 1;
    assert(std::abs(sin_value)<=1);
    double theta = 1/3.*std::asin(sin_value);
    double m = std::cos(theta) - 1./std::sqrt(3.)*std::sin(theta)*std::sin(phir);
    double beta = std::sin(phir)/m;
    double kappa = std::cos(phir)/m*c;
    double F = std::sqrt(s.J_2) + beta * s.I_1/3. - kappa;
    return F;
}


template <int DisplacementDim>
boost::optional<std::tuple<typename SCDamageModel2<DisplacementDim>::KelvinVector,
                           std::unique_ptr<typename MechanicsBase<
                               DisplacementDim>::MaterialStateVariables>,
                           typename SCDamageModel2<DisplacementDim>::KelvinMatrix>>
SCDamageModel2<DisplacementDim>::integrateStress(
    double const t, ProcessLib::SpatialPosition const& x, double const /*dt*/,
    KelvinVector const& eps_prev, KelvinVector const& eps,
    KelvinVector const& sigma_prev,
    typename MechanicsBase<DisplacementDim>::MaterialStateVariables const&
    /*material_state_variables*/,
    double const T) const
{
    const auto damaged = _damage_state(t, x)[0];
    const auto C = this->getElasticTensor(t, x, T, damaged);
    KelvinVector sigma_try = sigma_prev + C * (eps - eps_prev);

    return {std::make_tuple(sigma_try, createMaterialStateVariables(),
                            C)};
}


template <int DisplacementDim>
typename SCDamageModel2<DisplacementDim>::KelvinMatrix
SCDamageModel2<DisplacementDim>::getElasticTensor(
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

template class SCDamageModel2<2>;
template class SCDamageModel2<3>;

}
}  // namespace Solids
}  // namespace MaterialLib
