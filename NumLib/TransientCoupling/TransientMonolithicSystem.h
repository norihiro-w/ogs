/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 */

#pragma once

#include "NumLib/Coupling/MonolithicProblem.h"
#include "TransientCoupledSystem.h"

namespace NumLib
{

class AbstractTransientMonolithicSystem : public NumLib::AbstractMonolithicSystem<ITransientCoupledSystem>
{
public:
    virtual ~AbstractTransientMonolithicSystem() {}
};


}
