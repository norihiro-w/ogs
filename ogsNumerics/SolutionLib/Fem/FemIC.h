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

#include "BaseLib/CodingTools.h"
#include "GeoLib/GeoObject.h"
#include "MeshLib/Core/IMesh.h"
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
    FemIC(const MeshLib::IMesh* msh)
    : _msh(msh)
    {
    }

    ///
    virtual ~FemIC()
    {
        BaseLib::releaseObjectsInStdVector( _vec_func ); 
    }

    /// add a distribution
    void addDistribution(const GeoLib::GeoObject* geo, const NumLib::ITXFunction* ic_func);

    /// return the number of registered distributions
    size_t getNumberOfDistributions() const {return _vec_geo.size();};

    /// setup 
    void setup(NumLib::ITXDiscreteFunction<double> &u0) const;

private:
    const MeshLib::IMesh* _msh;
    std::vector<const GeoLib::GeoObject*> _vec_geo;
    std::vector<const NumLib::ITXFunction*> _vec_func;

};


}
