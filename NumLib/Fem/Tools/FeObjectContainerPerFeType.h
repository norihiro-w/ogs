/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 */

#pragma once

#include "BaseLib/CodingTools.h"

#include "NumLib/Fem/FiniteElement/IFemElement.h"
#include "IFeObjectContainer.h"
#include "FeObjectCachePerFeType.h"
#include "FemElementCatalog.h"
#include "IMeshElementToFemElementType.h"

namespace NumLib
{

/**
 * \brief Finite element object containers based on element shape types
 */
class FeObjectContainerPerFeType
 : public IFeObjectContainer
{
public:
    /**
     *
     * @param fe_catalog    FE type catalog
     * @param shape2feType  Link between shape type and FE type
     * @param msh           Mesh
     */
    FeObjectContainerPerFeType(const FemElementCatalog* fe_catalog, const IMeshElementToFemElementType* ele2feType, const MeshLib::Mesh* msh)
    : _cache(fe_catalog, msh), _ele2feType(ele2feType)
    {
    };

    /**
     * Copy constructor
     * @param src
     */
    FeObjectContainerPerFeType(const FeObjectContainerPerFeType &src)
    : _cache(src._cache), _ele2feType(src._ele2feType)
    {
    }

    /**
     *
     */
    virtual ~FeObjectContainerPerFeType() {};

    /**
     * return a copy
     * @return
     */
    virtual FeObjectContainerPerFeType* clone() const
    {
        return new FeObjectContainerPerFeType(*this);
    }

    /**
     *
     * @param e     Mesh element
     * @return a pointer to IFiniteElement object
     */
    virtual IFiniteElement* getFeObject(const MeshLib::Element &e) const
    {
        int fe_type = _ele2feType->getFeType(e, 1);//e.getCurrentOrder());
        IFiniteElement* fe = _cache.getFeObject(fe_type);
        fe->setMeshElement(e);
        return fe;
    }

private:
    const FeObjectCachePerFeType _cache;
    const IMeshElementToFemElementType* _ele2feType;
};



}
