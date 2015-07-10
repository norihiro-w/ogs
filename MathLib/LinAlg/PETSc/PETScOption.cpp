/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "PETScOption.h"

#include <petscksp.h>

namespace MathLib
{

PETScOption::PETScOption()
{
    solver_type = KSPBCGS;
    precon_type = PCNONE;
    max_iterations = 1e6;
    error_tolerance = 1.e-16;
}

PETScOption::SolverType PETScOption::getSolverType(const std::string &solver_name)
{
#define RETURN_TYPE(str, TypeName1, TypeName2) \
    if (#TypeName1==str) return TypeName2;

    RETURN_TYPE(solver_name, CG, KSPCG);
    RETURN_TYPE(solver_name, BiCGSTAB, KSPBCGS);
    RETURN_TYPE(solver_name, GMRES, KSPGMRES);

    return "";
#undef RETURN_SOLVER_ENUM_IF_SAME_STRING
}

PETScOption::PreconType PETScOption::getPreconType(const std::string &precon_name)
{
#define RETURN_TYPE(str, TypeName1, TypeName2) \
    if (#TypeName1==str) return TypeName2;

    RETURN_TYPE(precon_name, NONE, PCNONE);
    RETURN_TYPE(precon_name, JACOBI, PCJACOBI);
    RETURN_TYPE(precon_name, ILU, PCILU);

    return "";
#undef RETURN_SOLVER_ENUM_IF_SAME_STRING
}

} //MathLib
