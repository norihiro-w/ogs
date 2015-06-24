/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file FeObjectContainerPerElement.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <vector>

#include "BaseLib/CodingTools.h"
#include "MeshLib/Elements/Element.h"

#include "IFeObjectContainer.h"
#include "FemElementCatalog.h"

namespace NumLib
{
class IFiniteElement;

/**
 * \brief Finite element object containers
 */
class FeObjectContainerPerElement : public IFeObjectContainer
{
public:
    /**
     *
     * @param msh
     * @param ele_size
     */
    FeObjectContainerPerElement(const FemElementCatalog* fe_catalog, MeshLib::Mesh* msh)
    : _fe_catalog(fe_catalog), _msh(msh)
    {
        _vec_fem.resize(msh->getNElements());
    }

    /**
     * Copy constructor
     * @param src
     */
    FeObjectContainerPerElement(const FeObjectContainerPerElement &src)
    : _fe_catalog(src._fe_catalog), _msh(src._msh), _vec_fem(src._vec_fem)
    {

    }

    ///
    virtual ~FeObjectContainerPerElement()
    {
    	for (auto p : _vec_fem)
    		delete p;
    }

    /**
     *
     * @return
     */
    virtual FeObjectContainerPerElement* clone() const
    {
        return new FeObjectContainerPerElement(*this);
    }

    /**
     *
     * @param i
     * @param fe_type
     */
    void addFiniteElement(size_t i, int fe_type)
    {
        _vec_fem[i] = _fe_catalog->createFeObject(fe_type);
    }

    /**
     *
     * @param e
     * @return
     */
    virtual IFiniteElement* getFeObject(const MeshLib::Element &e) const
    {
        return _vec_fem[e.getID()];
    }

private:
    const FemElementCatalog* _fe_catalog;
    MeshLib::Mesh* _msh;
    std::vector<IFiniteElement*> _vec_fem;
};

}
