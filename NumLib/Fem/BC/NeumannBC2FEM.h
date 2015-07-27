/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 */

#pragma once

#include <vector>
#include <cstdlib>

namespace GeoLib
{
class GeoObject;
}

namespace MeshGeoToolsLib
{
class MeshNodeSearcher;
class BoundaryElementsSearcher;
}

namespace NumLib
{
class ITXFunction;
class IFeObjectContainer;

/**
 * Neumann BC
 */
class NeumannBC2FEM
{
public:
    /// 
    NeumannBC2FEM(MeshGeoToolsLib::MeshNodeSearcher& nodeSearcher, MeshGeoToolsLib::BoundaryElementsSearcher& beSearcher,
            const double &current_time,
            IFeObjectContainer &feObjects,
            const GeoLib::GeoObject &_geo, const NumLib::ITXFunction &_bc_func,
            std::vector<size_t> &_vec_nodes, std::vector<double> &_vec_values);
};



}
