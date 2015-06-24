/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 */

#pragma once

#include "BaseLib/CodingTools.h"

#include "NumLib/Fem/FiniteElement/IFemElement.h"
#include "FeObjectContainerPerFeType.h"
#include "FemElementCatalog.h"
#include "LagrangeFemElementCatalogBuilder.h"
#include "MeshElementShapeToFemElementType.h"

namespace NumLib
{

/**
 * \brief Lagrangian finite element object containers
 */
class LagrangeFeObjectContainer: public FeObjectContainerPerFeType
{
public:
    /**
     *
     * @param msh
     */
    explicit LagrangeFeObjectContainer(const MeshLib::Mesh* msh)
    : FeObjectContainerPerFeType(&_fe_catalog, &_shape2fetype, msh), _order(1)
    {
        LagrangeFemElementCatalogBuilder::construct(_fe_catalog, _shape2fetype);
    }

    /**
     * Copy constructor
     * @param src
     */
    LagrangeFeObjectContainer(const LagrangeFeObjectContainer &src)
    : FeObjectContainerPerFeType(src), _shape2fetype(src._shape2fetype), _order(src._order)
    {
    }

    /**
     *
     */
    virtual ~LagrangeFeObjectContainer() {}

    virtual LagrangeFeObjectContainer* clone() const
    {
        return new LagrangeFeObjectContainer(*this);
    }

    /**
     * set polynomial order
     * @param order
     */
    void setPolynomialOrder(size_t order) { _order = order; };

    /**
     * get FE object
     * @param e
     * @return
     */
    virtual IFiniteElement* getFeObject(const MeshLib::Element &e) const;

    /**
     * get FE object
     * @param e
     * @param order
     * @return
     */
    virtual IFiniteElement* getFeObject(const MeshLib::Element &e, size_t order) const;

private:
    FemElementCatalog _fe_catalog;
    MeshElementShapeToFemElementType _shape2fetype;
    size_t _order;
};



}
