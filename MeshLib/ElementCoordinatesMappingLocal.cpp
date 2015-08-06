/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "ElementCoordinatesMappingLocal.h"

#include <limits>
#include <cassert>

#include "GeoLib/AnalyticalGeometry.h"

#include "MeshLib/Elements/Element.h"
#include "MeshLib/Node.h"
#include "MathLib/MathTools.h"
#include "MathLib/Point3d.h"
#include "MathLib/Vector3.h"

namespace detail
{

/// rotate points to local coordinates
void rotateToLocal(
    const MathLib::RotationMatrix &matR2local,
    std::vector<MathLib::Point3d> &points)
{
    for (auto& p : points)
        p = matR2local*p;
}

/// get a rotation matrix to the global coordinates
/// it computes R in x=R*x' where x is original coordinates and x' is local coordinates
void getRotationMatrixToGlobal(
    const unsigned element_dimension,
    const unsigned global_dim,
    const std::vector<MathLib::Point3d> &points,
    MathLib::RotationMatrix &matR)
{
    // compute R in x=R*x' where x are original coordinates and x' are local coordinates
    if (element_dimension == 1) {
        MathLib::Vector3 xx(points[0], points[1]);
        xx.normalize();
        if (global_dim == 2)
            GeoLib::compute2DRotationMatrixToX(xx, matR);
        else
            GeoLib::compute3DRotationMatrixToX(xx, matR);
#ifdef OGS_USE_EIGEN
        matR.transposeInPlace();
#else
        matR.transpose();
#endif
    } else if (global_dim == 3 && element_dimension == 2) {
        // get plane normal
        MathLib::Vector3 plane_normal;
        double d;
        std::tie(plane_normal, d) = GeoLib::getNewellPlane(points);

        // compute a rotation matrix to XY
        GeoLib::computeRotationMatrixToXY(plane_normal, matR);
        // set a transposed matrix
#ifdef OGS_USE_EIGEN
        matR.transposeInPlace();
#else
        matR.transpose();
#endif
    }

}
}   // namespace detail

namespace MeshLib
{

ElementCoordinatesMappingLocal::ElementCoordinatesMappingLocal(
    const Element& e,
    const CoordinateSystem &global_coords)
: _coords(global_coords) //, _matR2global(3,3)
{
    assert(e.getDimension() <= global_coords.getDimension());
    _points.reserve(e.getNNodes());
    for(unsigned i = 0; i < e.getNNodes(); i++)
        _points.emplace_back(*(e.getNode(i)));

    auto const element_dimension = e.getDimension();
    auto const global_dimension = global_coords.getDimension();

    if (global_dimension == element_dimension)
    {
#if defined(OGS_USE_EIGEN)
        _matR2global.setIdentity();
#else
        _matR2global = 0;
        for (unsigned i=0; i<global_dimension; i++)
            _matR2global(i,i) = 1;
#endif
        return;
    }

    detail::getRotationMatrixToGlobal(element_dimension, global_dimension, _points, _matR2global);
#if defined(OGS_USE_EIGEN)
    detail::rotateToLocal(_matR2global.transpose(), _points);
#elif defined(OGS_USE_BLAZE)
    detail::rotateToLocal(blaze::trans(_matR2global), _points);
#else
    RotationMatrix* m(_matR2global.transpose());
    detail::rotateToLocal(*m, _points);
    delete m;
#endif
}

} // MeshLib
