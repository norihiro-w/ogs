/**
 * \file
 * \date   2015-01-16
 * \brief  Definition of the Point3d class.
 *
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef POINT3D_H_
#define POINT3D_H_

#include <limits>

#ifdef OGS_USE_EIGEN
#include <Eigen/Eigen>
#endif

#include "TemplatePoint.h"
#include "MathLib/LinAlg/Dense/DenseMatrix.h"

namespace MathLib
{
typedef MathLib::TemplatePoint<double,3> Point3d;
} // end namespace MathLib

bool operator< (MathLib::Point3d const & p0, MathLib::Point3d const & p1);

/**
 * lexicographic comparison of points
 */
bool operator<= (MathLib::Point3d const & p0, MathLib::Point3d const & p1);

/**
 * lexicographical comparison of points taking an epsilon into account
 * @param p0 first input Point3d
 * @param p1 second input Point3d
 * @param tol tolerance (if in the comparison operation the property fabs(p0[k] - p1[k]) < tol
 *     holds for the k-th coordinate the points are assumed the be equal in this coordinate)
 * @return true, if p0 is lexicographically smaller than p1
 */
bool lessEq(const MathLib::Point3d& p0,
            const MathLib::Point3d& p1,
            double tol = std::numeric_limits<double>::epsilon());

/**
 * rotation of points
 * @param mat a rotation matrix
 * @param p   a point to be transformed
 * @return a rotated point
 */
inline MathLib::Point3d operator*(const MathLib::DenseMatrix<double> &mat, const MathLib::Point3d &p)
{
    double* new_coords(mat*p.getCoords());
    MathLib::Point3d new_p(new_coords);
    delete [] new_coords;
    return new_p;
}

#ifdef OGS_USE_EIGEN
/**
 * rotation of points
 * @param mat a rotation matrix
 * @param p   a point to be transformed
 * @return a rotated point
 */
template <typename T_DERIVED>
MathLib::Point3d operator*(const Eigen::MatrixBase<T_DERIVED> &mat, const MathLib::Point3d &p)
{
    Eigen::Vector3d  x_new = mat * Eigen::Map<Eigen::Vector3d>(const_cast<double*>(p.getCoords()));
    return MathLib::Point3d(x_new.data());
}
#endif

#endif /* POINT3D_H_ */

