/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file Ogs6FemData.cpp
 *
 * Created on 2012-07-17 by Norihiro Watanabe
 */

#include "Ogs6FemData.h"

namespace THMCLib
{

Ogs6FemData* Ogs6FemData::_obj = 0;

Ogs6FemData* Ogs6FemData::getInstance()
{
    if (_obj==0) _obj = new Ogs6FemData();
    return _obj;
}

}
