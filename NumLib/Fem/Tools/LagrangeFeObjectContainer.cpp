/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file LagrangeFeObjectContainer.cpp
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#include "LagrangeFeObjectContainer.h"


namespace NumLib
{

IFiniteElement* LagrangeFeObjectContainer::getFeObject(const MeshLib::Element &e) const
{
    //e.setCurrentOrder(_order);
    IFiniteElement* fe = FeObjectContainerPerFeType::getFeObject(e);
    fe->setMeshElement(e);
    return fe;
}

IFiniteElement* LagrangeFeObjectContainer::getFeObject(const MeshLib::Element &e, size_t order) const
{
    //setPolynomialOrder(order);
    //e.setCurrentOrder(order);
    return getFeObject(e);
}


} //end
