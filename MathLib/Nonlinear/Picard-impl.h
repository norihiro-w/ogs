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


namespace MathLib
{

namespace Nonlinear
{

template<class T_FUNCTOR, class T_VALUE>
bool Picard::solve(T_FUNCTOR &functor, T_VALUE &x)
{
    T_VALUE x_old(x);
    T_VALUE dx(x);

    const bool checkAbsError = (_abs_tol<std::numeric_limits<double>::max());
    const bool checkRelError = (_rel_tol<std::numeric_limits<double>::max());

    INFO("------------------------------------------------------------------");
    INFO("*** PICARD nonlinear solver");
    INFO("-> iteration started");
    bool converged = false;
    std::size_t itr_cnt = 0;
    double x_norm = -1.;
    double abs_error = -1.;
    double rel_error = -1.;
    for (itr_cnt=0; itr_cnt<_max_itr; itr_cnt++) {
        functor(x_old, x);
        dx = x;
        dx -= x_old;

        abs_error = norm(dx, _normType);
        if (checkRelError) {
            x_norm = norm(x, _normType);
            if (x_norm>.0)
                rel_error = abs_error / x_norm;
            else
                rel_error = abs_error;
        }
        converged = (abs_error < _abs_tol && rel_error < _rel_tol);
        if (_printErrors)
            INFO("-> %d: ||dx||=%1.3e, ||x||=%1.3e, ||dx||/||x||=%1.3e", itr_cnt, abs_error, x_norm, rel_error);

#ifdef DEBUG_PICARD
        printout(std::cout, itr_cnt, x_new, dx);
#endif
        if (converged) {
            break;
        }
        x_old = x;
    }

    INFO("-> iteration finished");
    if (_max_itr==1) {
        INFO("status    : iteration not required");
    } else {
        INFO("status    : %s", (converged ? "CONVERGED" : "***DIVERGED***"));
    }
    INFO("iteration : %d/%d", itr_cnt, _max_itr);
    if (checkAbsError)
        INFO("abs error = %1.3e (tolerance=%1.3e)", abs_error, _abs_tol);
    if (checkRelError)
        INFO("rel error = %1.3e (tolerance=%1.3e)", rel_error, _rel_tol);
    INFO("norm type : %s", convertVecNormTypeToString(_normType).c_str());
    INFO("------------------------------------------------------------------");

    this->_n_iterations = itr_cnt;
    this->_abs_error = abs_error;
    this->_rel_error = rel_error;

    return converged;
}


template<class T_FUNCTOR>
bool Picard::solve(T_FUNCTOR &functor, MathLib::IVector &x)
{
    BaseLib::MPIEnvironment mpi;
    std::unique_ptr<MathLib::IVector> x_old(x.duplicate());
    std::unique_ptr<MathLib::IVector> dx(x.duplicate());

    const bool checkAbsError = (_abs_tol<std::numeric_limits<double>::max());
    const bool checkRelError = (_rel_tol<std::numeric_limits<double>::max());

    if (mpi.root()) {
        INFO("------------------------------------------------------------------");
        INFO("*** PICARD nonlinear solver");
        INFO("-> iteration started");
    }
    bool converged = false;
    std::size_t itr_cnt = 0;
    double x_norm = -1.;
    double abs_error = -1.;
    double rel_error = -1.;
    for (itr_cnt=0; itr_cnt<_max_itr; itr_cnt++) {
        functor(x);
        *dx = x;
        *dx -= *x_old;

        abs_error = norm(*dx, _normType);
        if (checkRelError) {
            x_norm = norm(x, _normType);
            if (x_norm>.0)
                rel_error = abs_error / x_norm;
            else
                rel_error = abs_error;
        }
        converged = (abs_error < _abs_tol && rel_error < _rel_tol);
        //if (_printErrors && mpi.root())
        INFO("-> %d: ||dx||=%1.3e, ||x||=%1.3e, ||dx||/||x||=%1.3e", itr_cnt, abs_error, x_norm, rel_error);

#ifdef DEBUG_PICARD
        printout(std::cout, itr_cnt, x_new, *dx);
#endif
        if (converged) {
            break;
        }
        *x_old = x;
    }

    if (mpi.root()) {
        INFO("-> iteration finished");
        if (_max_itr==1) {
            INFO("status    : iteration not required");
        } else {
            INFO("status    : %s", (converged ? "CONVERGED" : "***DIVERGED***"));
        }
        INFO("iteration : %d/%d", itr_cnt, _max_itr);
        if (checkAbsError)
            INFO("abs error = %1.3e (tolerance=%1.3e)", abs_error, _abs_tol);
        if (checkRelError)
            INFO("rel error = %1.3e (tolerance=%1.3e)", rel_error, _rel_tol);
        INFO("norm type : %s", convertVecNormTypeToString(_normType).c_str());
        INFO("------------------------------------------------------------------");
    }

    this->_n_iterations = itr_cnt;
    this->_abs_error = abs_error;
    this->_rel_error = rel_error;

    return converged;
}

#ifdef DEBUG_PICARD
template<class T_VALUE>
inline void Picard::printout(std::ostream& os, std::size_t i, T_VALUE& x_new, T_VALUE& dx)
{
    os << "-> " << i <<": x=(";
    for (std::size_t i=0; i<x_new.size(); i++)
        os << x_new[i] << " ";
    os << "), dx=(";
    for (std::size_t i=0; i<dx.size(); i++)
        os << dx[i] << " ";
    os << ")\n";
}

// in case of double
template<>
inline void Picard::printout(std::ostream& os, std::size_t i, double& x_new, double& dx)
{
    os << "-> " << i <<": x=" << x_new << ", dx=" << dx << "\n";
}
#endif

}

} //end
