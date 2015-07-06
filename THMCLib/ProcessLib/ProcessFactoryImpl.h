/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file ProcessFactoryImpl.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include "Process.h"
#include "ProcessFactoryBase.h"

namespace THMCLib
{

/**
 * \brief This class provides implementation of TeastFactoryBase interface.
 * It is used in TEST and TEST_F macros.
 */
template <class T>
class ProcessFactoryImpl : public ProcessFactoryBase
{
public:
  virtual Process* createProcess() { return new T; }
};

}
