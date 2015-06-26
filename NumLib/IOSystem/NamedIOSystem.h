/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file NamedIOSystem.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include "IOSystem.h"
#include "NamedIO.h"
#include "DynamicIOSystem.h"

namespace NumLib
{

/**
 *
 */
template <class T_BASE>
class NamedIOSystem : public NamedIO<IOSystem<T_BASE> >
{
public:
    virtual ~NamedIOSystem() {};
};

/**
 *
 */
template <class T_BASE>
class NamedDynamicIOSystem : public NamedIO<DynamicIOSystem<T_BASE> >
{
public:
    virtual ~NamedDynamicIOSystem() {};
};

} //end
