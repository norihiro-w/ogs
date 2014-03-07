/**
 * \file
 * \author Karsten Rink
 * \date   2013-07-05
 * \brief  Definition of mesh to geometry conversion.
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef CONVERTMESHTOGEO_H_
#define CONVERTMESHTOGEO_H_

#include <limits>

namespace GeoLib
{
class GEOObjects;
}

namespace MeshLib
{

	class Mesh;

	/**
	 * Converts a 2D mesh into a geometry.
	 * A new geometry with the name of the mesh will be inserted into geo_objects, consisting
	 * of points identical with mesh nodes and one surface representing the mesh. Triangles are
	 * converted to geometric triangles, quads are split into two triangles, all other elements
	 * are ignored.
	 */
	bool convertMeshToGeo(const MeshLib::Mesh &mesh, GeoLib::GEOObjects &geo_objects, double eps = std::numeric_limits<double>::epsilon());

} // namespace

#endif /* CONVERTMESHTOGEO_H_ */