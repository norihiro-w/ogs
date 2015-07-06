/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file IMedium.h
 *
 * Created on 2012-11-29 by Norihiro Watanabe
 */

#pragma once

#include "MediumType.h"

namespace MaterialLib
{

/**
 * Interface to any kinds of medium
 */
struct IMedium
{
    virtual ~IMedium() {};
    virtual MediumType getMediumType() const = 0;
};

} //end
 
