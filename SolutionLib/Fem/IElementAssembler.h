/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 */

#pragma once

#include "MathLib/DataType.h"

namespace MeshLib
{
class Element;
};

namespace AssemblerLib
{
class LocalToGlobalIndexMap;
}

namespace NumLib
{
class TimeStep;
}

namespace SolutionLib
{

/**
 * \brief Interface class for element local assembler for transient case
 */
class IElementAssembler
{
public:
    virtual ~IElementAssembler() {}

    virtual void reset(const MeshLib::Element &e) = 0;
    //virtual void reset(const MeshLib::Element &e, const AssemblerLib::LocalToGlobalIndexMap &localDof) = 0;

    virtual void linear(  const NumLib::TimeStep &time,
                            const MathLib::LocalVector &local_u_n1, const MathLib::LocalVector &local_u_n, 
                            MathLib::LocalMatrix &localA, MathLib::LocalVector &localRhs) = 0;

    virtual void residual(  const NumLib::TimeStep &time,
                            const MathLib::LocalVector &local_u_n1, const MathLib::LocalVector &local_u_n,
                            MathLib::LocalVector &localR) = 0;

    virtual void jacobian(  const NumLib::TimeStep &time,
                            const MathLib::LocalVector &local_u_n1, const MathLib::LocalVector &local_u_n,
                            MathLib::LocalMatrix &localJ) = 0;
};

} //end
