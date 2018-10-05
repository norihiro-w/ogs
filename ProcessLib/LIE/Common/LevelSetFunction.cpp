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

template <class T_VEC, class T_V>
inline bool inList(T_VEC const& vec, T_V const& v)
{
	return (std::find(vec.begin(), vec.end(), v) != vec.end());
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
	std::unordered_map<int,int> const& fracID_to_local,
	Eigen::Vector3d const& x)
{
	//pre-calculate levelsets for all fractures
	std::vector<double> levelsets(frac_props.size());
	for (std::size_t i=0; i<frac_props.size(); i++)
		levelsets[i] = Heaviside(levelset_fracture(*frac_props[i], x));

	std::vector<double> enrichments(frac_props.size() + junction_props.size());
	//fractures possibly with branches
	for (std::size_t i=0; i<frac_props.size(); i++)
	{
		auto const* frac = frac_props[i];
		double enrich = levelsets[i];
		for (std::size_t j=0; j<frac->branches.size(); j++)
			enrich *= Heaviside(levelset_branch(*frac->branches[j], x));
		enrichments[i] = enrich;
	}

	// junctions
	for (std::size_t i=0; i<junction_props.size(); i++)
	{
		auto const* junction = junction_props[i];
		auto fid1 = fracID_to_local.at(junction->fracture_IDs[0]);
		auto fid2 = fracID_to_local.at(junction->fracture_IDs[1]);
		double enrich = levelsets[fid1]*levelsets[fid2];
		enrichments[i + frac_props.size()] *= enrich;
	}

	return enrichments;
}


std::vector<double> du_global_enrichments(
	std::size_t this_frac_index,
	std::vector<FractureProperty*> const& frac_props,
	std::vector<JunctionProperty*> const& junction_props,
	std::unordered_map<int,int> const& fracID_to_local,
	Eigen::Vector3d const& x)
{
	auto const& this_frac = *frac_props[this_frac_index];
	//pre-calculate levelsets for all fractures
	std::vector<double> levelsets(frac_props.size());
	for (std::size_t i=0; i<frac_props.size(); i++)
		levelsets[i] = Heaviside(levelset_fracture(*frac_props[i], x));

	std::vector<double> enrichments(frac_props.size() + junction_props.size());
	enrichments[this_frac_index] = 1.0;

	//fractures possibly with branches
	if (frac_props.size()>1)
	{
		for (auto const& branch : this_frac.branches)
		{
			if (branch->master_fracture_ID != this_frac.fracture_id)
				continue;

			double singned = boost::math::sign(this_frac.normal_vector.dot(branch->normal_vector_branch));
			auto slave_fid = fracID_to_local.at(branch->slave_fracture_ID);
			double enrich = singned * Heaviside(levelsets[slave_fid]);
			enrichments[slave_fid] = enrich;
		}
	}

	// junctions
	for (unsigned i=0; i<junction_props.size(); i++)
	{
		auto const* junction = junction_props[i];
		if (!inList(junction->fracture_IDs, this_frac.fracture_id))
			continue;

		auto another_frac_id = (junction->fracture_IDs[0]==this_frac.fracture_id) ? junction->fracture_IDs[1] : junction->fracture_IDs[0];
		auto fid = fracID_to_local.at(another_frac_id);
		double enrich = Heaviside(levelsets[fid]);
		enrichments[i + frac_props.size()] *= enrich;
	}


	return enrichments;
}

}  // namespace LIE
}  // namespace ProcessLib
