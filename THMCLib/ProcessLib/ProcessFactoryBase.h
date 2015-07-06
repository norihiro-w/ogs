/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file ProcessFactoryBase.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include "BaseLib/CodingTools.h"
#include "Process.h"


namespace THMCLib
{

/**
 * \brief Defines the abstract factory interface that creates instances
 * of a Test object.
 */
class ProcessFactoryBase
{
public:
  virtual ~ProcessFactoryBase() {}

  /// Creates a test instance to run. The instance is both created and destroyed
  /// within TestInfoImpl::Run()
  virtual Process* createProcess() = 0;

protected:
  ProcessFactoryBase() {}

private:
  DISALLOW_COPY_AND_ASSIGN(ProcessFactoryBase);
};

}
