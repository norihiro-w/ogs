/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file IC2FEM.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
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
 * \brief A class constructing I.C. data for FEM
 */
class IC2FEM
{
public:
    ///
    IC2FEM(const MeshLib::Mesh &msh, const GeoLib::GeoObject &geo, const NumLib::ITXFunction &ic_func, std::vector<size_t> &vec_nodes, std::vector<double> &vec_values);
};


}
