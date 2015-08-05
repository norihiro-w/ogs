/*!
   \file  PETScTools.h
   \brief Declaration of a function related to PETSc solver interface to assign
         the Dirichlet boundary conditions.

   \author Wenqing Wang
   \version
   \date Nov 2011 - Sep 2013

   \copyright
    Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
               Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license
*/

#ifndef PETSCTOOLS_H_
#define PETSCTOOLS_H_

#include <vector>

#include "PETScMatrix.h"
#include "PETScVector.h"

namespace MathLib
{
/*!
   \brief apply known solutions to a system of linear equations

   \param A                 Coefficient matrix
   \param b                 RHS vector
   \param vec_knownX_id    a vector of known solution entry IDs
   \param vec_knownX_x     a vector of known solutions
   \param penalty_scaling value for scaling some matrix and right hand side
        entries to enforce some conditions
*/
void applyKnownSolution(PETScMatrix &A, PETScVector &b,
                        const std::vector<std::size_t> &vec_knownX_id,
                        const std::vector<double> &vec_knownX_x);
} // end of namespace MathLib

#endif //end  of PETSCTOOLS_H_

