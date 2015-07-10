/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef ILINEARSOLVER_H_
#define ILINEARSOLVER_H_

#include <boost/property_tree/ptree.hpp>

#include "LinAlgLibType.h"

namespace MathLib
{

class IVector;
class IMatrix;

class ILinearSolver
{
public:
    virtual ~ILinearSolver() {}

    /// returns a library type
    virtual LinAlgLibType getLinAlgLibType() const = 0;

    /**
     * configure linear solvers
     * @param option
     */
    virtual void setOption(const boost::property_tree::ptree &option) = 0;

    /**
     * solve a given linear equations
     *
     * @param b     RHS vector
     * @param x     Solution vector
     */
    virtual void solve(MathLib::IVector &b, MathLib::IVector &x) = 0;

    /// apply prescribed values to a system of linear equations
    virtual void imposeKnownSolution(IMatrix &A, IVector &b, const std::vector<std::size_t> &vec_knownX_id,
    		const std::vector<double> &vec_knownX_x, double penalty_scaling = 1e+10) = 0;
};

} // MathLib

#endif //ILINEARSOLVER_H_

