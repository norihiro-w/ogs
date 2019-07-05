/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <cmath>
#include <memory>
#include <tuple>

#include "MathLib/KelvinVector.h"
#include "LinearElasticIsotropic.h"
#include "ProcessLib/Parameter/Parameter.h"

namespace MaterialLib
{
namespace Solids
{

namespace SCDamage
{

/**
 * \brief A class for computing the SCDamageModel model
 *
 */
template <int DisplacementDim>
class SCDamageModel final : public MechanicsBase<DisplacementDim>
{
public:
    struct MaterialProperties
    {
        using P = ProcessLib::Parameter<double>;

        MaterialProperties(
            P const& youngs_modulus, P const& youngs_modulus_damaged, P const& poissons_ratio,
            P const& friction_angle, P const& cohesion)
            : _youngs_modulus(youngs_modulus), _youngs_modulus_damaged(youngs_modulus_damaged), _poissons_ratio(poissons_ratio),
              _friction_angle(friction_angle), _cohesion(cohesion)
        {
        }

        P const& _youngs_modulus;
        P const& _youngs_modulus_damaged;
        P const& _poissons_ratio;
        P const& _friction_angle;
        P const& _cohesion;
    };


    struct MaterialStateVariables
        : public MechanicsBase<DisplacementDim>::MaterialStateVariables
    {
        void pushBackState() override {}
        MaterialStateVariables& operator=(MaterialStateVariables const&) =
            default;
        typename MechanicsBase<DisplacementDim>::MaterialStateVariables&
        operator=(typename MechanicsBase<DisplacementDim>::
                      MaterialStateVariables const& state) noexcept override
        {
            return operator=(static_cast<MaterialStateVariables const&>(state));
        }
    };

    static int const KelvinVectorSize =
        MathLib::KelvinVector::KelvinVectorDimensions<DisplacementDim>::value;
    using ResidualVectorType = Eigen::Matrix<double, KelvinVectorSize, 1>;
    using JacobianMatrix = Eigen::Matrix<double, KelvinVectorSize,
                                         KelvinVectorSize, Eigen::RowMajor>;

    using KelvinVector =
        MathLib::KelvinVector::KelvinVectorType<DisplacementDim>;
    using KelvinMatrix =
        MathLib::KelvinVector::KelvinMatrixType<DisplacementDim>;

    using Parameter = ProcessLib::Parameter<double>;

    std::unique_ptr<
        typename MechanicsBase<DisplacementDim>::MaterialStateVariables>
    createMaterialStateVariables() const override
    {
        return std::unique_ptr<
            typename MechanicsBase<DisplacementDim>::MaterialStateVariables>{
            new MaterialStateVariables};
    }

    SCDamageModel(
        typename SCDamageModel<DisplacementDim>::MaterialProperties mp,
        Parameter const& damage_state, bool const check_MC)
        :  _mp(std::move(mp)),
          _damage_state(damage_state),
          _check_MC(check_MC)
    {
    }

    double computeFreeEnergyDensity(
        double const /*t*/,
        ProcessLib::SpatialPosition const& /*x*/,
        double const /*dt*/,
        KelvinVector const& eps,
        KelvinVector const& sigma,
        typename MechanicsBase<DisplacementDim>::
            MaterialStateVariables const& /* material_state_variables */)
        const override
    {
        return eps.dot(sigma) / 2;
    }

    double yieldFunctionMC(
        double const /*t*/,
        ProcessLib::SpatialPosition const& /*x*/,
        double const /*dt*/,
        KelvinVector const& sigma,
        typename MechanicsBase<DisplacementDim>::MaterialStateVariables const&
            material_state_variables) const;

    boost::optional<std::tuple<KelvinVector,
                               std::unique_ptr<typename MechanicsBase<
                                   DisplacementDim>::MaterialStateVariables>,
                               KelvinMatrix>>
    integrateStress(
        double const t, ProcessLib::SpatialPosition const& x, double const dt,
        KelvinVector const& eps_prev, KelvinVector const& eps,
        KelvinVector const& sigma_prev,
        typename MechanicsBase<DisplacementDim>::MaterialStateVariables const&
            material_state_variables,
        double const T) const override;

    ConstitutiveModel getConstitutiveModel() const override
    {
        return ConstitutiveModel::SCDamageModel;
    }

    KelvinMatrix getElasticTensor(double const t,
                                  ProcessLib::SpatialPosition const& x,
                                  double const T, double const damaged) const;

    MaterialProperties getMaterialProperties() const { return _mp; }

private:
    MaterialProperties _mp;
    Parameter const& _damage_state;
    bool const _check_MC;
};

extern template class SCDamageModel<2>;
extern template class SCDamageModel<3>;

}
}  // end of namespace Solids
}  // namespace MaterialLib
