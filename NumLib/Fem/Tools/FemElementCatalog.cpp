/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file FemElementFactory.cpp
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#include "FemElementCatalog.h"


namespace NumLib
{

IFiniteElement* FemElementCatalog::createFeObject(const int fe_id) const
{
    if (!exist(fe_id))
        return NULL;

    std::map<int, FemElementFactoryBase*>::const_iterator itr = _map_fe_factories.find(fe_id);
    return itr->second->create();
}


}
