/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PETSC_OPTION_H_
#define PETSC_OPTION_H_

#include <string>

namespace MathLib
{

/**
 * \brief Option for PETSc solver
 */
struct PETScOption
{
    typedef std::string SolverType;
    typedef std::string PreconType;

    /// Linear solver type
    SolverType solver_type;
    /// Preconditioner type
    PreconType precon_type;
    /// Maximum iteration count
    long max_iterations;
    /// Error tolerance
    double error_tolerance;
    /// Extra option
    std::string extra_arg;
    /// Arguments for solver and preconditioner. This variable is always preferred
    /// to other variables.
    std::string solver_precon_arg;

    /**
     * Constructor
     *
     * Default options are CG, no preconditioner, iteration count 500 and
     * tolerance 1e-10. Default matrix storage type is CRS.
     */
    PETScOption();

    /// Destructor
    ~PETScOption() {}

    /**
     * return a linear solver type from the solver name
     *
     * @param solver_name
     * @return a linear solver type
     *      If there is no solver type matched with the given name, INVALID
     *      is returned.
     */
    static SolverType getSolverType(const std::string &solver_name);

    /**
     * return a preconditioner type from the name
     *
     * @param precon_name
     * @return a preconditioner type
     *      If there is no preconditioner type matched with the given name, NONE
     *      is returned.
     */
    static PreconType getPreconType(const std::string &precon_name);
};

}
#endif //PETSC_OPTION_H_
