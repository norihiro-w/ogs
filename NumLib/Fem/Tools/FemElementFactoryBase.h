/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file ProcessFactoryBase.h
 *
 * Created on 2012-11-27 by Norihiro Watanabe
 */

#pragma once

#include "BaseLib/CodingTools.h"
#include "NumLib/Fem/FiniteElement/IFemElement.h"

namespace NumLib
{

/**
 * \brief Defines the abstract factory interface that creates instances
 * of a Test object.
 */
class FemElementFactoryBase
{
public:
  virtual ~FemElementFactoryBase() {}

  /// Creates a test instance to run. The instance is both created and destroyed
  /// within TestInfoImpl::Run()
  virtual IFiniteElement* create() = 0;

protected:
  FemElementFactoryBase() {}

private:
  DISALLOW_COPY_AND_ASSIGN(FemElementFactoryBase);
};

}
