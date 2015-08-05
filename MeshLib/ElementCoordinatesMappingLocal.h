/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#ifndef ELEMENTCOORDINATESMAPPINGLOCAL_H_
#define ELEMENTCOORDINATESMAPPINGLOCAL_H_

#include <vector>

#ifdef OGS_USE_EIGEN
#include <Eigen/Eigen>
#else
#include "MathLib/LinAlg/Dense/DenseMatrix.h"
#endif

#include "MathLib/Point3d.h"
#include "MathLib/DataType.h"

#include "MeshLib/CoordinateSystem.h"

namespace MeshLib
{
    class Element;
}

namespace MeshLib
{

/**
 * This class maps node coordinates on intrinsic coordinates of the given element.
 */
class ElementCoordinatesMappingLocal final
{
public:
    /**
     * Constructor
     * \param e                     Mesh element whose node coordinates are mapped
     * \param global_coord_system   Global coordinate system
     */
    ElementCoordinatesMappingLocal(const Element &e, const CoordinateSystem &global_coord_system);

    /// return the global coordinate system
    const CoordinateSystem getGlobalCoordinateSystem() const { return _coords; }

    /// return mapped coordinates of the node
    MathLib::Point3d const& getMappedCoordinates(std::size_t node_id) const
    {
        return _points[node_id];
    }

    /// return a rotation matrix converting to global coordinates
    const MathLib::RotationMatrix& getRotationMatrixToGlobal() const {return _matR2global;}

private:
    const CoordinateSystem _coords;
    std::vector<MathLib::Point3d> _points;
    MathLib::RotationMatrix _matR2global;
};

}

#endif // ELEMENTCOORDINATESMAPPINGLOCAL_H_
