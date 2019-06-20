/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <vector>

#include "NumLib/Extrapolation/ExtrapolatableElement.h"
#include "ProcessLib/LocalAssemblerInterface.h"
#include "ProcessLib/Parameter/Parameter.h"

namespace ProcessLib
{
namespace SmallDeformationWithPT
{
struct SmallDeformationWithPTLocalAssemblerInterface
    : public ProcessLib::LocalAssemblerInterface,
      public NumLib::ExtrapolatableElement
{
    virtual std::size_t setIPDataInitialConditions(
        std::string const& name, double const* values,
        int const integration_order) = 0;

    virtual std::size_t setIPDataInitialConditions(
        std::string const& name, Parameter<double> const& function) = 0;

    virtual std::vector<double> getSigma() const = 0;

    virtual std::vector<double> getEpsilon() const = 0;

    virtual std::vector<double> getEpsilonMechanical() const = 0;

    virtual std::vector<double> const& getIntPtSigma(
        const double /*t*/,
        GlobalVector const& /*current_solution*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> const& getIntPtEpsilon(
        const double /*t*/,
        GlobalVector const& /*current_solution*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& cache) const = 0;

    virtual MeshLib::Element const& getMeshElement() const = 0;

    virtual unsigned getNumberOfIntegrationPoints() const = 0;

    virtual void assembleResidual(
        double const t,
        std::vector<double>& local_rhs_data) const = 0;
};

}  // namespace SmallDeformationWithPT
}  // namespace ProcessLib
