/*!
   \file  PETScTools.cpp
   \brief Definition of a function related to PETSc solver interface to assign
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

#include "PETScTools.h"

namespace MathLib
{

void applyKnownSolution(PETScMatrix &A, PETScVector &b,
                        const std::vector<std::size_t>  &vec_knownX_id,
                        const std::vector<double> &vec_knownX_x)
{
    A.finalizeAssembly();

    A.setRowsColumnsZero(vec_knownX_id);
    A.finalizeAssembly();

    b.finalizeAssembly();
    if(vec_knownX_id.size() > 0)
    {
        b.set(vec_knownX_id, vec_knownX_x);
    }

    b.finalizeAssembly();
}

} // end of namespace MathLib


