/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "LevelSetFunction.h"

#include <boost/math/special_functions/sign.hpp>

#include "BranchProperty.h"
#include "FractureProperty.h"
#include "JunctionProperty.h"

namespace
{
// Heaviside step function
inline double Heaviside(double v)
{
    return (v < 0.0) ? 0.0 : 1.0;
}

}  // namespace

namespace ProcessLib
{
namespace LIE
{

double levelset_fracture(FractureProperty const& frac, Eigen::Vector3d const& x)
{
    return boost::math::sign(frac.normal_vector.dot(x - frac.point_on_fracture));
}

double levelset_branch(BranchProperty const& branch, Eigen::Vector3d const& x)
{
	return boost::math::sign(branch.normal_vector_branch.dot(x - branch.coords));
}


std::vector<double> u_global_enrichments(
	std::vector<FractureProperty*> const& frac_props,
	std::vector<JunctionProperty*> const& junction_props,
	Eigen::Vector3d const& x)
{
	//pre-calculate levelsets for all fractures
	std::vector<double> levelsets(frac_props.size());
	for (auto const* frac : frac_props)
		levelsets[frac->fracture_id] = Heaviside(levelset_fracture(*frac, x));

	std::vector<double> enrichments(frac_props.size() + junction_props.size());
	//fractures possibly with branches
	for (auto const* frac : frac_props)
	{
		double enrich = levelsets[frac->fracture_id];
		for (std::size_t i=0; i<frac->branches.size(); i++)
			enrich *= Heaviside(levelset_branch(*frac->branches[i], x));
		enrichments[frac->fracture_id] = enrich;
	}

	// junctions
	for (auto const* junction : junction_props)
	{
		double enrich = levelsets[junction->fracture_IDs[0]]*levelsets[junction->fracture_IDs[1]];
		enrichments[junction->junction_id + frac_props.size()] *= enrich;
	}

	return enrichments;
}


std::vector<double> du_global_enrichments(
	FractureProperty const& this_frac,
	std::vector<FractureProperty*> const& frac_props,
	std::vector<JunctionProperty*> const& junction_props,
	Eigen::Vector3d const& x)
{

	//pre-calculate levelsets for all fractures
	std::vector<double> levelsets(frac_props.size());
	for (auto const* frac : frac_props)
		levelsets[frac->fracture_id] = Heaviside(levelset_fracture(*frac, x));

	std::vector<double> enrichments(frac_props.size() + junction_props.size());
	enrichments[this_frac.fracture_id] = 1.0;

	//fractures possibly with branches
	for (auto const* frac : frac_props)
	{
		for (auto const* branch : frac->branches)
		{
			if (branch->master_fracture_ID != this_frac.fracture_id)
				continue;

			double singned = boost::math::sign(this_frac.normal_vector.dot(branch->normal_vector_branch));
			double enrich = singned * Heaviside(levelsets[branch->slave_fracture_ID]);
			enrichments[branch->slave_fracture_ID] = enrich;
		}
	}

	// junctions
	for (auto const* junction : junction_props)
	{
		if (std::find(junction->fracture_IDs.begin(), junction->fracture_IDs.end(), this_frac.fracture_id)
			== junction->fracture_IDs.end())
			continue;

		auto other_frac_id = (junction->fracture_IDs[0]==this_frac.fracture_id) ? junction->fracture_IDs[1] : junction->fracture_IDs[0];
		double enrich = Heaviside(levelsets[other_frac_id]);
		enrichments[junction->junction_id + frac_props.size()] *= enrich;
	}


	return enrichments;
}

}  // namespace LIE
}  // namespace ProcessLib
