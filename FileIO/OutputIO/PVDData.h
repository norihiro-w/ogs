/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file PVDData.h
 *
 * Created on 2012-07-17 by Norihiro Watanabe
 */

#pragma once

#include <string>
#include <vector>


struct PVDData
{
    struct VTK_Info
    {
        double timestep;
        std::string vtk_file;
    };

    std::vector<VTK_Info> vec_dataset;
};
