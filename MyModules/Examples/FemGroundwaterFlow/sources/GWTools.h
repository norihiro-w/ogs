/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef GWTOOLS_H
#define GWTOOLS_H

#include <logog/include/logog.hpp>

#include "NumLib/Function/FemNodalFunction.h"
#include "NumLib/Function/FemIntegrationPointFunction.h"

void getDarcyVelocity(const NumLib::FemNodalFunctionScalar &f_p, NumLib::FEMIntegrationPointFunctionVector &vel);

#endif
