/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef EIGENLISLINEARSOLVER_H_
#define EIGENLISLINEARSOLVER_H_

#include <vector>

#include <boost/property_tree/ptree_fwd.hpp>
#include <lis.h>

#include "MathLib/LinAlg/ILinearSolver.h"
#include "MathLib/LinAlg/Lis/LisOption.h"

namespace MathLib
{
class EigenVector;
class EigenMatrix;

/**
 * Linear solver using Lis library with Eigen matrix and vector objects
 */
class EigenLisLinearSolver final : public ILinearSolver
{
public:
    /**
     * Constructor
     * @param A         Coefficient matrix object
     * @param option    A pointer to a linear solver option. In case you omit
     *                  this argument, default settings follow those of
     *                  LisOption struct.
     */
    EigenLisLinearSolver(EigenMatrix &A, boost::property_tree::ptree const*const option = nullptr);

    LinAlgLibType getLinAlgLibType() const {return LinAlgLibType::EigenLis;}

    /**
     * parse linear solvers configuration
     */
    void setOption(const boost::property_tree::ptree &option);

    /**
     * copy linear solvers options
     */
    void setOption(const LisOption &option) { _option = option; }

    /**
     * get linear solver options
     */
    LisOption &getOption() { return _option; }

    /**
     * solve a given linear equations
     *
     * @param b     RHS vector
     * @param x     Solution vector
     */
    void solve(IVector &b, IVector &x);

    /// apply prescribed values to a system of linear equations
    void imposeKnownSolution(IMatrix &A, IVector &b, const std::vector<std::size_t> &vec_knownX_id,
            const std::vector<double> &vec_knownX_x, double penalty_scaling = 1e+10);

private:
    EigenMatrix& _A;
    LisOption _option;
};

} // MathLib

#endif //EIGENLISLINEARSOLVER_H_

