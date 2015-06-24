/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 */

#include "FeObjectCachePerFeType.h"

namespace NumLib
{

FeObjectCachePerFeType::FeObjectCachePerFeType(const FeObjectCachePerFeType& src)
: _fe_catalog(src._fe_catalog), _msh(src._msh)
{
    std::map<int, IFiniteElement*>::const_iterator itr = src._mapFeObj.begin();
    for (; itr!=src._mapFeObj.end(); ++itr) {
        _mapFeObj.insert(*itr);
    }
}

IFiniteElement* FeObjectCachePerFeType::getFeObject(const int fe_type) const
{
    IFiniteElement *fe = 0;
    if (_mapFeObj.count(fe_type)==0) {
        fe = _fe_catalog->createFeObject(fe_type);
        _mapFeObj[fe_type] = fe;
    } else {
        fe = _mapFeObj[fe_type];
    }
    return fe;
}

}
