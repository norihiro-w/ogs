/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file FeObjectCachePerFeType.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <map>

#include "BaseLib/CodingTools.h"
#include "MeshLib/Mesh.h"

#include "NumLib/Fem/FiniteElement/IFemElement.h"
#include "FemElementCatalog.h"

namespace NumLib
{

/**
 * \brief Cache system for finite element objects 
 */
class FeObjectCachePerFeType
{
public:
    /**
     *
     * @param fe_catalog
     * @param msh
     */
    explicit FeObjectCachePerFeType(const FemElementCatalog* fe_catalog, const MeshLib::Mesh* msh)
    : _fe_catalog(fe_catalog), _msh(msh) {};

    /**
     *
     * @param src
     */
    FeObjectCachePerFeType(const FeObjectCachePerFeType& src);

    /**
     *
     */
    virtual ~FeObjectCachePerFeType()
    {
        for (auto itr=_mapFeObj.begin(); itr!=_mapFeObj.end(); ++itr)
            if (itr->second!=0) delete itr->second;
    }

    /**
     * get a finite element object
     * @param fe_type
     * @return
     */
    IFiniteElement* getFeObject(const int fe_type) const;

    /**
     * get a pointer to this mesh
     * @return
     */
    //const MeshLib::Mesh* getMesh() const {return _msh;};

    /**
     * get a pointer to FE catalog
     * @return
     */
    const FemElementCatalog* getCatalog() const {return _fe_catalog;};

private:
    const FemElementCatalog* _fe_catalog;
    const MeshLib::Mesh* _msh;
    mutable std::map<int, IFiniteElement*> _mapFeObj;
};

}
