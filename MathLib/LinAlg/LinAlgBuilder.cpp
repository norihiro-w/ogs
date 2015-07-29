/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "LinAlgBuilder.h"

#include "IVector.h"
#include "IMatrix.h"
#include "ILinearSolver.h"

#include "Dense/DenseVector.h"
#include "Dense/GlobalDenseMatrix.h"
#include "Eigen/EigenMatrix.h"
#include "Eigen/EigenVector.h"
#include "Eigen/EigenLinearSolver.h"
#ifdef USE_LIS
#include "Lis/LisVector.h"
#include "Lis/LisMatrix.h"
#include "Lis/LisLinearSolver.h"
#include "EigenLis/EigenLisLinearSolver.h"
#endif
#ifdef USE_PETSC
#include "PETSc/PETScVector.h"
#include "PETSc/PETScMatrix.h"
#include "PETSc/PETScLinearSolver.h"
#endif

namespace MathLib
{

IVector* LinAlgBuilder::duplicateVector(IVector &v)
{
    switch (v.getLinAlgLibType())
    {
    case LinAlgLibType::Dense:
        return new MathLib::DenseVector<double>(static_cast<MathLib::DenseVector<double>&>(v));
    case LinAlgLibType::Eigen:
        return new MathLib::EigenVector(static_cast<MathLib::EigenVector&>(v));
#ifdef USE_LIS
    case LinAlgLibType::EigenLis:
        return new MathLib::EigenVector(static_cast<MathLib::EigenVector&>(v));
    case LinAlgLibType::Lis:
        return new MathLib::LisVector(static_cast<MathLib::LisVector&>(v));
#endif
#ifdef USE_PETSC
    case LinAlgLibType::PETSc:
        return new MathLib::PETScVector(static_cast<MathLib::PETScVector&>(v));
#endif
    default:
        return nullptr;
    }
}

IVector* LinAlgBuilder::generateVector(LinAlgLibType libType, std::size_t n, bool is_global_size, std::vector<std::size_t> const* ghost_ids)
{
    switch (libType)
    {
    case LinAlgLibType::Dense:
        return new MathLib::DenseVector<double>(n);
    case LinAlgLibType::Eigen:
        return new MathLib::EigenVector(n);
#ifdef USE_LIS
    case LinAlgLibType::EigenLis:
        return new MathLib::EigenVector(n);
    case LinAlgLibType::Lis:
        return new MathLib::LisVector(n);
#endif
#ifdef USE_PETSC
    case LinAlgLibType::PETSc:
        return new MathLib::PETScVector(n, is_global_size, ghost_ids);
#endif
    default:
        return nullptr;
    }
}

IMatrix* LinAlgBuilder::generateMatrix(LinAlgLibType libType, std::size_t n, const MatrixOption* opt)
{
    switch (libType)
    {
//	case LinAlgLibType::Dense:
//		return new MathLib::GlobalDenseMatrix<double>(n);
    case LinAlgLibType::Eigen:
        return new MathLib::EigenMatrix(n);
#ifdef USE_LIS
    case LinAlgLibType::EigenLis:
        return new MathLib::EigenMatrix(n);
    case LinAlgLibType::Lis:
        return new MathLib::LisMatrix(n);
#endif
#ifdef USE_PETSC
    case LinAlgLibType::PETSc:
        return opt ? new MathLib::PETScMatrix(n, *opt) : new MathLib::PETScMatrix(n);
#endif
    default:
        return nullptr;
    }

}

ILinearSolver* LinAlgBuilder::generateLinearSolver(LinAlgLibType libType, IMatrix* A, boost::property_tree::ptree const*const option)
{
    switch (libType)
    {
//	case LinAlgLibType::Dense:
//		return new MathLib::GlobalDenseMatrix<double>(n);
    case LinAlgLibType::Eigen:
        return new MathLib::EigenLinearSolver(*static_cast<EigenMatrix*>(A), option);
#ifdef USE_LIS
    case LinAlgLibType::EigenLis:
        return new MathLib::EigenLisLinearSolver(*static_cast<EigenMatrix*>(A), option);
    case LinAlgLibType::Lis:
        return new MathLib::LisLinearSolver(*static_cast<LisMatrix*>(A), option);
#endif
#ifdef USE_PETSC
    case LinAlgLibType::PETSc:
        return new MathLib::PETScLinearSolver(*static_cast<MathLib::PETScMatrix*>(A), option);
#endif
    default:
        return nullptr;
    }

}

} // MathLib


