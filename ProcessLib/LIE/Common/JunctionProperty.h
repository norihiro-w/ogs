/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <Eigen/Eigen>

namespace ProcessLib
{
namespace LIE
{

struct JunctionProperty
{
	int junction_id;
	int node_id;
	Eigen::Vector3d coords;
	std::array<int, 2> fracture_IDs;

    virtual ~JunctionProperty() = default;
};


}  // namespace LIE
}  // namespace ProcessLib
