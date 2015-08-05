/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "EigenLinearSolver.h"

#include <boost/property_tree/ptree.hpp>
#include <logog/include/logog.hpp>

#include "EigenVector.h"
#include "EigenMatrix.h"
#include "EigenTools.h"

namespace MathLib
{

using boost::property_tree::ptree;

namespace details
{

/// Template class for Eigen direct linear solvers
template <class T_SOLVER, class T_BASE>
class EigenDirectLinearSolver final : public T_BASE
{
public:
    explicit EigenDirectLinearSolver(EigenMatrix::RawMatrixType &A) : _A(A)
    {
        INFO("-> initialize with the coefficient matrix");
        _solver.compute(A);
        if(_solver.info()!=Eigen::Success) {
            ERR("Failed during Eigen linear solver initialization");
            return;
        }
    }

    void solve(EigenVector::RawVectorType &b, EigenVector::RawVectorType &x, EigenOption &/*opt*/) override
    {
        INFO("-> solve");
        x = _solver.solve(b);
        if(_solver.info()!=Eigen::Success) {
            ERR("Failed during Eigen linear solve");
            return;
        }
    }

private:
    T_SOLVER _solver;
    EigenMatrix::RawMatrixType& _A;
};

/// Template class for Eigen iterative linear solvers
template <class T_SOLVER, class T_BASE>
class EigenIterativeLinearSolver final : public T_BASE
{
public:
    explicit EigenIterativeLinearSolver(EigenMatrix::RawMatrixType &A) : _A(A)
    {
        INFO("-> initialize with the coefficient matrix");
        _solver.compute(A);
        if(_solver.info()!=Eigen::Success) {
            ERR("Failed during Eigen linear solver initialization");
            return;
        }
    }

    void solve(EigenVector::RawVectorType &b, EigenVector::RawVectorType &x, EigenOption &opt) override
    {
        INFO("-> solve");
        _solver.setTolerance(opt.error_tolerance);
        _solver.setMaxIterations(opt.max_iterations);
        x = _solver.solveWithGuess(b, x);
        if(_solver.info()!=Eigen::Success) {
            ERR("Failed during Eigen linear solve");
            return;
        }
        INFO("\t iteration: %d/%ld", _solver.iterations(), opt.max_iterations);
        INFO("\t residual: %e\n", _solver.error());
    }

private:
    T_SOLVER _solver;
    EigenMatrix::RawMatrixType& _A;
};

} // details

EigenLinearSolver::EigenLinearSolver(EigenMatrix &A, ptree const*const option)
{
    if (option)
        setOption(*option);

    A.getRawMatrix().makeCompressed();
    if (_option.solver_type==EigenOption::SolverType::SparseLU) {
        using SolverType = Eigen::SparseLU<EigenMatrix::RawMatrixType, Eigen::COLAMDOrdering<int>>;
        _solver = new details::EigenDirectLinearSolver<SolverType, IEigenSolver>(A.getRawMatrix());
    } else if (_option.solver_type==EigenOption::SolverType::BiCGSTAB) {
        using SolverType = Eigen::BiCGSTAB<EigenMatrix::RawMatrixType, Eigen::DiagonalPreconditioner<double>>;
        _solver = new details::EigenIterativeLinearSolver<SolverType, IEigenSolver>(A.getRawMatrix());
    } else if (_option.solver_type==EigenOption::SolverType::CG) {
        using SolverType = Eigen::ConjugateGradient<EigenMatrix::RawMatrixType, Eigen::Lower, Eigen::DiagonalPreconditioner<double>>;
        _solver = new details::EigenIterativeLinearSolver<SolverType, IEigenSolver>(A.getRawMatrix());
    }
}

void EigenLinearSolver::setOption(const ptree &option)
{
    boost::optional<ptree> ptSolver = option.get_child("LinearSolver");
    if (!ptSolver)
        return;

    boost::optional<std::string> solver_type = ptSolver->get_optional<std::string>("solver_type");
    if (solver_type) {
        _option.solver_type = _option.getSolverType(*solver_type);
    }
    boost::optional<std::string> precon_type = ptSolver->get_optional<std::string>("precon_type");
    if (precon_type) {
        _option.precon_type = _option.getPreconType(*precon_type);
    }
    boost::optional<double> error_tolerance = ptSolver->get_optional<double>("error_tolerance");
    if (error_tolerance) {
        _option.error_tolerance = *error_tolerance;
    }
    boost::optional<int> max_iteration_step = ptSolver->get_optional<int>("max_iteration_step");
    if (max_iteration_step) {
        _option.max_iterations = *max_iteration_step;
    }
}

void EigenLinearSolver::solve(IVector &b_, IVector &x_)
{
    auto &b = static_cast<EigenVector&>(b_);
    auto &x = static_cast<EigenVector&>(x_);

    INFO("------------------------------------------------------------------");
    INFO("*** Eigen solver computation");
    _solver->solve(b.getRawVector(), x.getRawVector(), _option);
    INFO("------------------------------------------------------------------");
}

void EigenLinearSolver::imposeKnownSolution(IMatrix &A, IVector &b, const std::vector<std::size_t> &vec_knownX_id,
		const std::vector<double> &vec_knownX_x, double penalty_scaling)
{
	applyKnownSolution(static_cast<EigenMatrix&>(A), static_cast<EigenVector&>(b), vec_knownX_id, vec_knownX_x, penalty_scaling);
}

} //MathLib
