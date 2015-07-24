/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef LINALGBUILDER_H_
#define LINALGBUILDER_H_

#include <cstdlib>
#include <vector>

#include <boost/property_tree/ptree.hpp>

#include "LinAlgLibType.h"
#include "MatrixOption.h"

namespace MathLib
{

class IVector;
class IMatrix;
class ILinearSolver;

class LinAlgBuilder
{
public:
    static IVector* duplicateVector(IVector &v);
    static IVector* generateVector(LinAlgLibType libType, std::size_t n, bool is_global_size = true, std::vector<std::size_t> const* ghost_ids = nullptr);
    static IMatrix* generateMatrix(LinAlgLibType libType, std::size_t n, const MatrixOption* opt = nullptr);
    static ILinearSolver* generateLinearSolver(LinAlgLibType libType, IMatrix*A, boost::property_tree::ptree const*const option = nullptr);
};

} // MathLib

#endif //LINALGBUILDER_H_

