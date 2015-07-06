/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file FemIC.h
 *
 * Created on 2012-09-22 by Norihiro Watanabe
 */

#pragma once

#include <vector>

#include "GeoLib/GeoObject.h"
#include "MeshLib/Mesh.h"
#include "NumLib/Function/ITXFunction.h"
#include "NumLib/Function/ITXDiscreteFunction.h"

namespace SolutionLib
{

/**
 * \brief IC data for FEM
 *
 * - mesh
 * - geometry
 * - distribution
 */
class FemIC
{
public:
    /**
     *
     * @param msh
     */
    FemIC(const MeshLib::Mesh* msh)
    : _msh(msh)
    {
    }

    ///
    virtual ~FemIC()
    {
    	for (auto p :_vec_func) delete p;
    }

    /// add a distribution
    void addDistribution(const GeoLib::GeoObject* geo, const NumLib::ITXFunction* ic_func);

    /// return the number of registered distributions
    size_t getNumberOfDistributions() const {return _vec_geo.size();};

    /// setup 
    void setup(MathLib::IVector &u0) const;

private:
    const MeshLib::Mesh* _msh;
    std::vector<const GeoLib::GeoObject*> _vec_geo;
    std::vector<const NumLib::ITXFunction*> _vec_func;

};


}
