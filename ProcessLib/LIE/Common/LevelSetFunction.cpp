/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "LevelSetFunction.h"

#include <boost/math/special_functions/sign.hpp>

#include "FractureProperty.h"

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
double calculateLevelSetFunction(FractureProperty const& frac, double const* x_)
{
    Eigen::Map<Eigen::Vector3d const> x(x_, 3);
    return Heaviside(
        boost::math::sign(frac.normal_vector.dot(x - frac.point_on_fracture)));
}

typedef struct BranchInfo
{
	std::vector<int> list_connected_discontinuity;
	long node_id;
	Eigen::Vector3d branch_vector;
} t_BranchInfo;

std::vector<long> vec_junction_nodes;
std::vector<std::pair<int, int> > vec_junction_discon_id;
std::vector<t_BranchInfo> vec_branches;

bool isConnectingToThisBranch(int this_discon_id, t_BranchInfo &this_branch)
{
	bool is_connected_to_this_branch = false;
	for (int k = 0; k<this_branch.list_connected_discontinuity.size(); k++)
	{
		if (this_branch.list_connected_discontinuity[k] == this_discon_id)
			is_connected_to_this_branch = true;
	}
	return is_connected_to_this_branch;
}

double levelSet4Branch(FractureProperty const& frac, 
					   std::vector<double> const& local_levelsets, double const* x_)
{
	Eigen::Map<Eigen::Vector3d const> x(x_, 3);

	// for branches: psi_b_i(x) = sign(n_a*n_ab_i)*H(f_b_i)
	double global_levelset = local_levelsets[frac.fracture_id];
	for (t_BranchInfo &branch : vec_branches)
	{
		if (branch.list_connected_discontinuity[1] != frac.fracture_id)
			continue; // should be branch

		// Get the signed distance from the branch point to the reference point
		Eigen::Vector3d point_on_branch; // = (double*)m_msh->nod_vector[this_branch.node_id]->GetCoordinates();
		Eigen::Vector3d rel_vector = x - point_on_branch;
		double singned = boost::math::sign(branch.branch_vector.dot(rel_vector));
		// Update the level set function
		global_levelset *= Heaviside(singned);
	}

	return global_levelset;
}


void levelSet4Junction(std::vector<double> const& local_levelsets, double const* x_)
{
	std::vector<double> vec_global_levelset_fct;

	Eigen::Map<Eigen::Vector3d const> x(x_, 3);

	int this_discon_id;
	FractureProperty* frac;

	int nr_discontinuity_junctions;
	int nr_discontinuous_features;
	// for junctions: %psi_r_i = H(f_r_i)
	for (int i=0; i<nr_discontinuity_junctions; i++) 
	{
		const int junction_delta_id = nr_discontinuous_features + i;
		// Junction information
		int dis1_id = vec_junction_discon_id[i].first;
		int dis2_id = vec_junction_discon_id[i].second;
		if (dis1_id != this_discon_id && dis2_id != this_discon_id)
			continue;

		int conneced_discon_id = dis1_id;
		if (dis1_id == this_discon_id)
			conneced_discon_id = dis2_id;

		// Set the level set
		vec_global_levelset_fct[junction_delta_id] = Heaviside(local_levelsets[conneced_discon_id]);
	}
}

}  // namespace LIE
}  // namespace ProcessLib
