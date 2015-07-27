/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <limits>
#include <memory>

#include <logog/include/logog.hpp>

#include "BaseLib/MPITools.h"
#include "MathLib/LinAlg/IVector.h"
#include "MathLib/LinAlg/IMatrix.h"
#include "MathLib/LinAlg/ILinearSolver.h"
#include "MathLib/LinAlg/LinAlgBuilder.h"

namespace MathLib
{

namespace Nonlinear
{

template<class F_RESIDUAL, class F_DX, class T_VALUE>
bool NewtonRaphson::solve(F_RESIDUAL &f_residual, F_DX &f_dx, T_VALUE &x)
{
    const bool checkAbsResidualError = (_r_abs_tol < std::numeric_limits<double>::max());
    const bool checkRelResidualError = (_r_rel_tol < std::numeric_limits<double>::max());
    const bool checkRelDxError = (_dx_rel_tol > .0);
    const bool needXNorm =  (checkRelResidualError || checkRelDxError);

    BaseLib::MPIEnvironment mpi;
    if (mpi.root()) {
        INFO("------------------------------------------------------------------");
        INFO("*** NEWTON-RAPHSON nonlinear solver");
        INFO("-> iteration started");
    }

    // evaluate initial residual
    T_VALUE r(x);
    f_residual(x, r);

    // check convergence
    double r_norm = norm(r, _normType);

    T_VALUE dx(x);
    double dx_norm = norm(dx, _normType);

    double x_norm = -1.;
    if (needXNorm)
        x_norm = norm(x, _normType);

    bool converged = ((r_norm < _r_abs_tol && r_norm < _r_rel_tol*x_norm)
                     || (checkRelDxError && dx_norm < _dx_rel_tol*x_norm));
    if (_printErrors && mpi.root())
        INFO("-> %d: ||r||=%1.3e, ||dx||=%1.3e, ||x||=%1.3e, ||dx||/||x||=%1.3e", 0, r_norm, dx_norm, x_norm, x_norm==0 ? dx_norm : dx_norm/x_norm);

    std::size_t itr_cnt = 0;
    if (!converged) {
        for (itr_cnt=1; itr_cnt<_max_itr; itr_cnt++) {
            // solve dx=-J^-1*r
            f_dx(x, r, dx);
            x += dx;
            // evaluate residual
            f_residual(x, r);
#ifdef DEBUG_NEWTON_RAPHSON
            printout(std::cout, itr_cnt, x, r, dx);
#endif
            // check convergence
            r_norm = norm(r, _normType);
            dx_norm = norm(dx, _normType);
            if (needXNorm)
                x_norm = norm(x, _normType);
            converged = ((r_norm < _r_abs_tol && r_norm < _r_rel_tol*x_norm)
                        || (checkRelDxError && dx_norm < _dx_rel_tol*x_norm));
            if (_printErrors && mpi.root())
                INFO("-> %d: ||r||=%1.3e, ||dx||=%1.3e, ||x||=%1.3e, ||dx||/||x||=%1.3e", itr_cnt, r_norm, dx_norm, x_norm, x_norm==0 ? dx_norm : dx_norm/x_norm);

            if (converged)
                break;
        }
    }

    if (mpi.root()) {
        INFO("-> iteration finished");
        if (_max_itr==1) {
            INFO("status    : iteration not required");
        } else {
            INFO("status    : %s", (converged ? "CONVERGED" : "***DIVERGED***"));
        }
        INFO("iteration : %d/%d", itr_cnt, _max_itr);
        if (checkAbsResidualError)
            INFO("abs. res. : %1.3e (tolerance=%1.3e)", r_norm, _r_abs_tol);
        if (checkRelResidualError)
            INFO("rel. res. : %1.3e (tolerance=%1.3e)", x_norm==0?r_norm:r_norm/x_norm, _r_rel_tol);
        if (checkRelDxError)
            INFO("dx        : %1.3e (tolerance=%1.3e)", x_norm==0?dx_norm:dx_norm/x_norm, _dx_rel_tol);
        INFO("norm type : %s", convertVecNormTypeToString(_normType).c_str());
        INFO("------------------------------------------------------------------");
    }

    this->_n_iterations = itr_cnt;
    this->_r_abs_error = r_norm;
    this->_r_rel_error = (x_norm==0 ? r_norm : r_norm / x_norm);
    this->_dx_rel_error = (x_norm==0 ? dx_norm : dx_norm / x_norm);

    return converged;
}

template<class F_RESIDUAL, class F_J>
bool NewtonRaphson::solve(F_RESIDUAL &f_residual, F_J &f_jacobian, IVector &x)
{
    const bool checkAbsResidualError = (_r_abs_tol < std::numeric_limits<double>::max());
    const bool checkRelResidualError = (_r_rel_tol < std::numeric_limits<double>::max());
    const bool checkRelDxError = (_dx_rel_tol > .0);
    const bool needXNorm =  (checkRelResidualError || checkRelDxError);

    BaseLib::MPIEnvironment mpi;
    if (mpi.root()) {
        INFO("------------------------------------------------------------------");
        INFO("*** NEWTON-RAPHSON nonlinear solver");
        INFO("-> iteration started");
    }

    // evaluate initial residual
    if (_r == nullptr)
        _r = x.duplicate();
    IVector &r(*_r);
    f_residual(x, r);

    // check convergence
    double r_norm = norm(r, _normType);

    std::unique_ptr<IVector> dxx(x.duplicate());
    IVector &dx(*dxx);
    double dx_norm = norm(dx, _normType);

    double x_norm = -1.;
    if (needXNorm)
        x_norm = norm(x, _normType);

    bool converged = ((r_norm < _r_abs_tol && r_norm < _r_rel_tol*x_norm)
                     || (checkRelDxError && dx_norm < _dx_rel_tol*x_norm));
    if (_printErrors && mpi.root())
        INFO("-> %d: ||r||=%1.3e, ||dx||=%1.3e, ||x||=%1.3e, ||dx||/||x||=%1.3e", 0, r_norm, dx_norm, x_norm, x_norm==0 ? dx_norm : dx_norm/x_norm);

    std::size_t itr_cnt = 0;
    if (!converged) {
        for (itr_cnt=1; itr_cnt<(_max_itr+1); itr_cnt++) {
            // form J
            f_jacobian(x, *_J);
            // solve dx=-J^-1*r
            *_b = r;
            *_b *= -1;
            auto ls = LinAlgBuilder::generateLinearSolver(x.getLinAlgLibType(), _J, _option);
            ls->solve(*_b, dx);
            delete ls;

            x += dx;
            x.assemble();
            //x.write("x" + std::to_string(itr_cnt) + ".txt");
            // evaluate residual
            f_residual(x, r);
#ifdef DEBUG_NEWTON_RAPHSON
            mpi.barrier();
            printout(std::cout, itr_cnt, x, r, dx);
#endif
            // check convergence
            r_norm = norm(r, _normType);
            dx_norm = norm(dx, _normType);
            if (needXNorm)
                x_norm = norm(x, _normType);
            converged = ((r_norm < _r_abs_tol && r_norm < _r_rel_tol*x_norm)
                        || (checkRelDxError && dx_norm < _dx_rel_tol*x_norm));
            if (_printErrors && mpi.root())
                INFO("-> %d: ||r||=%1.3e, ||dx||=%1.3e, ||x||=%1.3e, ||dx||/||x||=%1.3e", itr_cnt, r_norm, dx_norm, x_norm, x_norm==0 ? dx_norm : dx_norm/x_norm);

            if (converged)
                break;
        }
    }

    if (mpi.root()) {
        INFO("-> iteration finished");
        if (_max_itr==1) {
            INFO("status    : iteration not required");
        } else {
            INFO("status    : %s", (converged ? "CONVERGED" : "***DIVERGED***"));
        }
        INFO("iteration : %d/%d", itr_cnt, _max_itr);
        if (checkAbsResidualError)
            INFO("abs. res. : %1.3e (tolerance=%1.3e)", r_norm, _r_abs_tol);
        if (checkRelResidualError)
            INFO("rel. res. : %1.3e (tolerance=%1.3e)", x_norm==0?r_norm:r_norm/x_norm, _r_rel_tol);
        if (checkRelDxError)
            INFO("dx        : %1.3e (tolerance=%1.3e)", x_norm==0?dx_norm:dx_norm/x_norm, _dx_rel_tol);
        INFO("norm type : %s", convertVecNormTypeToString(_normType).c_str());
        INFO("------------------------------------------------------------------");
    }

    this->_n_iterations = itr_cnt;
    this->_r_abs_error = r_norm;
    this->_r_rel_error = (x_norm==0 ? r_norm : r_norm / x_norm);
    this->_dx_rel_error = (x_norm==0 ? dx_norm : dx_norm / x_norm);

    return converged;
}

#ifdef DEBUG_NEWTON_RAPHSON
template<class T_VALUE>
inline void NewtonRaphson::printout(std::ostream& os, std::size_t i, T_VALUE& x_new, T_VALUE& r, T_VALUE& dx)
{
    os << "-> " << i <<": x=(";
    for (std::size_t i=0; i<x_new.size(); i++)
        os << x_new[i] << " ";
    os << "), r=(";
    for (std::size_t i=0; i<dx.size(); i++)
        os << r[i] << " ";
    os << "), dx=(";
    for (std::size_t i=0; i<dx.size(); i++)
        os << dx[i] << " ";
    os << ")\n";
}

template<>
inline void NewtonRaphson::printout(std::ostream& os, std::size_t i, IVector& x_new, IVector& r, IVector& dx)
{
    os << "-> " << i <<": x=(";
    for (std::size_t i=x_new.getRangeBegin(); i<x_new.getRangeEnd(); i++)
        os << x_new[i] << " ";
    os << "), r=(";
    for (std::size_t i=r.getRangeBegin(); i<r.getRangeEnd(); i++)
        os << r[i] << " ";
    os << "), dx=(";
    for (std::size_t i=dx.getRangeBegin(); i<dx.getRangeEnd(); i++)
        os << dx[i] << " ";
    os << ")\n";
}

// in case of double
template<>
inline void NewtonRaphson::printout(std::ostream& os, std::size_t i, double& x_new, double& r, double& dx)
{
    os << "-> " << i <<": x=" << x_new << ", r=" << r << ", dx=" << dx << "\n";
}
#endif

}

} //end
