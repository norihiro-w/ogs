/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file TransientCouplingStructureBuilder.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include "NumLib/Coupling/CouplingStrucutreBuilder.h"
#include "NumLib/Coupling/Algorithm/TransientPartitionedAlgorithmFactory.h"

#include "TransientCoupledSystem.h"
#include "TransientMonolithicSystem.h"
#include "AsyncPartitionedSystem.h"


namespace NumLib
{

typedef class TemplateCouplingStrucutreBuilder
    <
    ITransientCoupledSystem,
    AbstractTransientMonolithicSystem,
    AsyncPartitionedSystem,
    TransientPartitionedAlgorithmFactory
    > TransientCoulplingStrucutreBuilder;

} //end
