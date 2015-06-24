/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file LagrangeFemElementCatalogBuilder.h
 *
 * Created on 2012-11-27 by Norihiro Watanabe
 */

#pragma once

namespace NumLib
{
class FemElementCatalog;
class MeshElementShapeToFemElementType;

/**
 * setup Lagrange FE type catalog
 */
class LagrangeFemElementCatalogBuilder
{
public:
    static void construct(FemElementCatalog &feCatalog, MeshElementShapeToFemElementType &shape2fetype);
};

}
