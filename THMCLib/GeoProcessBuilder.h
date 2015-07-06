/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file GeoProcessBuilder.h
 *
 * Created on 2012-07-13 by Norihiro Watanabe
 */

#pragma once

#include "ProcessLib/AbstractProcessBuilder.h"

/**
 * \brief Process builder for OGS6
 *
 * This class follows singleton pattern.
 */
class GeoProcessBuilder: public THMCLib::AbstractProcessBuilder
{
public:
    static THMCLib::AbstractProcessBuilder* getInstance();
private:
    static THMCLib::AbstractProcessBuilder* _obj;

public:
    virtual ~GeoProcessBuilder() {};

private:
    GeoProcessBuilder();
};
