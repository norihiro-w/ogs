/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file FemElementFactoryImpl.h
 *
 * Created on 2012-11-27 by Norihiro Watanabe
 */

#pragma once

#include "FemElementFactoryBase.h"

namespace NumLib
{

/**
 * \brief This class provides implementation of TeastFactoryBase interface.
 * It is used in TEST and TEST_F macros.
 */
template <class T>
class FemElementFactoryImpl : public FemElementFactoryBase
{
public:
  virtual IFiniteElement* create() { return new T(); }
};

}
