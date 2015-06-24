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

namespace MeshLib
{
class Mesh;
}

namespace NumLib
{
class ITXFunction;

/**
 * \brief A class constructing DirichletBC data for FEM
 */
class DirichletBC2FEM
{
public:
    ///
    DirichletBC2FEM(const MeshLib::Mesh &msh, const GeoLib::GeoObject &geo, const NumLib::ITXFunction &bc_func, std::vector<size_t> &vec_nodes, std::vector<double> &vec_values);
};


}
