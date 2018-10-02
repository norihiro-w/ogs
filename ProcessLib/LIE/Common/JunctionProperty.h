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
	std::size_t junction_id;
	std::size_t node_id;
	Eigen::Vector3d coords;
	std::size_t fracture1_ID;
	std::size_t fracture2_ID;

    virtual ~JunctionProperty() = default;
};


}  // namespace LIE
}  // namespace ProcessLib
