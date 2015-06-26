/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file ITXDiscreteFunction.h
 *
 * Created on 2012-08-17 by Norihiro Watanabe
 */

#pragma once

#include "ITXFunction.h"
#include "MathLib/LinAlg/LinAlgEnums.h"

namespace MathLib
{
class IVector;
}

namespace NumLib
{

/**
 * \brief Interface of any functions in space-time domain
 *
 * This class aims to be an abstract of spatially and temporally distributed data such as
 * physical quantity (e.g. head, stress) and material property (e.g. permeability).
 *
 * TXFunction
 * - is evaluated at particular position in space-time domain
 * - returns scalar or vector value
 * - has some attributes (e.g constant)
 *
 */
class ITXDiscreteFunction : public ITXFunction
{
public:
    virtual ~ITXDiscreteFunction() {}

    virtual double diff_norm(const ITXDiscreteFunction &v, MathLib::VecNormType normType) const = 0;
};

} //end


