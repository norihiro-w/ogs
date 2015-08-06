/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */


#pragma once

#ifdef OGS_USE_BLAZE

#include <blaze/Math.h>

namespace MathLib
{
/// Local dense matrix type (row-majored)
typedef blaze::DynamicMatrix<double,blaze::rowMajor> LocalMatrix;
/// Local dense vector type
typedef blaze::DynamicVector<double> LocalVector;

typedef blaze::StaticMatrix<double, 3u, 3u, blaze::rowMajor> RotationMatrix;

}

#else

#include <Eigen/Eigen>

namespace MathLib
{
/// Local dense matrix type (row-majored)
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> LocalMatrix;
/// Local dense vector type
typedef Eigen::VectorXd LocalVector;

#ifdef OGS_USE_EIGEN
typedef Eigen::Matrix<double, 3u, 3u, Eigen::RowMajor> RotationMatrix;
#else
typedef MathLib::DenseMatrix<double> RotationMatrix;
#endif

}

#endif
