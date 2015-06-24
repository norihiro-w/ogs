/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file Ogs5FemIO.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <string>
#include <vector>

#include "Ogs5FemData.h"

namespace ogs5
{

namespace Ogs5FemIO
{
bool read(const std::string &proj_path, Ogs5FemData &ogs5data);
};

}
