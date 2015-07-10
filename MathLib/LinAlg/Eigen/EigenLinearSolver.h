/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef EIGENLINEARSOLVER_H_
#define EIGENLINEARSOLVER_H_

#include <vector>

#include <boost/property_tree/ptree_fwd.hpp>

#include "MathLib/LinAlg/ILinearSolver.h"
#include "EigenVector.h"
#include "EigenOption.h"

namespace MathLib
{

class EigenMatrix;

class EigenLinearSolver final :public ILinearSolver
{
public:
    EigenLinearSolver(EigenMatrix &A, boost::property_tree::ptree const*const option = nullptr);

    virtual ~EigenLinearSolver()
    {
        delete _solver;
    }

    LinAlgLibType getLinAlgLibType() const {return LinAlgLibType::Eigen;}

    /**
     * parse linear solvers configuration
     */
    void setOption(const boost::property_tree::ptree &option);

    /**
     * copy linear solvers options
     */
    void setOption(const EigenOption &option) { _option = option; }

    /**
     * get linear solver options
     */
    EigenOption &getOption() { return _option; }

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

protected:
    class IEigenSolver
    {
    public:
        virtual ~IEigenSolver() = default;
        /**
         * execute a linear solver
         */
        virtual void solve(EigenVector::RawVectorType &b, EigenVector::RawVectorType &x, EigenOption &) = 0;
    };

    EigenOption _option;
    IEigenSolver* _solver;
};

} // MathLib

#endif //EIGENLINEARSOLVER_H_

