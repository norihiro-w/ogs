/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <vector>

#include <Eigen/Eigen>

namespace ProcessLib
{
namespace LIE
{
struct FractureProperty;
struct JunctionProperty;

// /// calculate the level set function
// /// \f$ \psi(\mathbf{x}) = H(|\mathbf{x} - \mathbf{x_d}|
// /// \mathrm{sign}[\mathbf{n_d} \cdot (\mathbf{x}-\mathbf{x_d}]) \f$ where
// /// \f$H(u)\f$ is the Heaviside step function, \f$\mathbf{x_d}\f$ is a point on
// /// the fracture plane, and \f$\mathbf{n_d}\f$ is the normal vector of a
// /// fracture plane
// double calculateLevelSetFunction(FractureProperty const& fracture_property,
//                                  double const* x);

double levelset_fracture(FractureProperty const& frac, Eigen::Vector3d const& x);

std::vector<double> u_global_enrichments(
	std::vector<FractureProperty*> const& frac_props,
	std::vector<JunctionProperty*> const& junction_props,
	Eigen::Vector3d const& x);

std::vector<double> du_global_enrichments(
	std::size_t this_frac_index,
	std::vector<FractureProperty*> const& frac_props,
	std::vector<JunctionProperty*> const& junction_props,
	Eigen::Vector3d const& x);

}  // namespace LIE
}  // namespace ProcessLib
