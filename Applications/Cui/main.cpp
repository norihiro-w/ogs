/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file main.cpp
 *
 * Created on 2012-07-17 by Norihiro Watanabe
 */

#include <iostream>
#include <exception>
#include <logog/include/logog.hpp>

#include "THMCLib/THMCSimulator.h"

int main ( int argc, char *argv[] )
{
    ogs6::ogsInit(argc, argv);
    auto sim = new ogs6::THMCSimulator(argc, argv);
    int returncode = sim->execute();
	delete sim;
    ogs6::ogsExit();

    return returncode;
}
