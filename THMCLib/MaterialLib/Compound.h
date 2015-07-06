/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file Compound.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <string>


namespace MaterialLib
{

typedef TemplateBasicModel<double> MolecularDiffusion;

struct Compound
{
    std::string name;
    bool is_mobile;
    MolecularDiffusion* molecular_diffusion;

    Compound()
    : is_mobile(true), molecular_diffusion(nullptr)
    {
    }

    ~Compound()
    {
    	delete molecular_diffusion;
    }

};

} //end

