/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file OutputBuilder.cpp
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#include "OutputBuilder.h"

#include "PVDOutput.h"
#include "TecplotOutput.h"

IOutput* OutputBuilder::create(const std::string &name)
{
    if (name.compare("PVD")==0) {
        return new PVDOutput();
    } else if (name.compare("TECPLOT")==0) {
        return new TecplotOutput();
    }
    return NULL;
}
