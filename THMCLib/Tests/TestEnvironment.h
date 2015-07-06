/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file TestEnvironment.h
 *
 * Created on 2012-12-10 by Norihiro Watanabe
 */

#pragma once

#include <gtest/gtest.h>

class TestEnvironment : public ::testing::Environment
{
 public:
  virtual ~TestEnvironment() {}
  // Override this to define how to set up the environment.
  virtual void SetUp()
  {
      //TODO delete old output data
  }
  // Override this to define how to tear down the environment.
  virtual void TearDown() {}
};

