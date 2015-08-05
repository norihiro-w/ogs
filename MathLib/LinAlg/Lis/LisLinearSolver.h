/**
 * \file
 * \author Norihiro Watanabe
 * \date   2013-04-16
 * \brief  Definition of the LisLinearSolver class.
 *
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef LISLINEARSOLVER_H_
#define LISLINEARSOLVER_H_

#include <vector>
#include <boost/property_tree/ptree.hpp>
#include <lis.h>

#include "MathLib/LinAlg/ILinearSolver.h"
#include "LisOption.h"
#include "LisVector.h"
#include "LisMatrix.h"

namespace MathLib
{

/**
 * \brief Linear solver using Lis (http://www.ssisc.org/lis/)
 *
 */
class LisLinearSolver :public ILinearSolver
{
public:
    /**
     * Constructor
     * @param A         Coefficient matrix object
     * @param option    A pointer to a linear solver option. In case you omit
     *                  this argument, default settings follow those of
     *                  LisOption struct.
     */
    LisLinearSolver(LisMatrix &A, boost::property_tree::ptree const*const option = nullptr);

    virtual ~LisLinearSolver() {}

    LinAlgLibType getLinAlgLibType() const {return LinAlgLibType::Lis;}

    /**
     * configure linear solvers
     * @param option
     */
    void setOption(const boost::property_tree::ptree &option);

    /**
     * configure linear solvers
     * @param option
     */
    void setOption(const LisOption &option) { _option = option; }

    /**
     * get linear solver options
     * @return
     */
    LisOption &getOption() { return _option; }

    /**
     * solve a given linear equations
     *
     * @param b     RHS vector
     * @param x     Solution vector
     */
    void solve(MathLib::IVector &b, MathLib::IVector &x);

    /// apply prescribed values to a system of linear equations
    void imposeKnownSolution(IMatrix &A, IVector &b, const std::vector<std::size_t> &vec_knownX_id,
    		const std::vector<double> &vec_knownX_x, double penalty_scaling = 1e+10);

private:
    LisMatrix& _A;
    LisOption _option;
};

} // MathLib

#endif //LISLINEARSOLVER_H_

