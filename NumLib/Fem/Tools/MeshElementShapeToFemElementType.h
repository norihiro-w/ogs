/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file MeshElementShapeToFemElementType.h
 *
 * Created on 2012-11-27 by Norihiro Watanabe
 */

#pragma once

#include <map>
#include "MeshLib/MeshEnums.h"
#include "MeshLib/Elements/Element.h"
#include "IMeshElementToFemElementType.h"

namespace NumLib
{

/**
 * \brief Converter from element shape to FE type
 */
class MeshElementShapeToFemElementType
: public IMeshElementToFemElementType
{
    typedef std::pair<MeshLib::MeshElemType, int> MyKey;
    typedef std::map<MyKey, int> MyTable;
public:
    virtual ~MeshElementShapeToFemElementType() {}

    /**
     * get FE type
     * @param e         Mesh element
     * @param order     Polynomial order
     * @return FE type
     */
    virtual int getFeType(const MeshLib::Element& e, const int order) const
    {
        MyKey key(e.getGeomType(), order);
        MyTable::const_iterator itr = _map_shape2fetype.find(key);
        if (itr != _map_shape2fetype.end())
            return itr->second;
        else
            return -1;
    }

    /**
     * add a mapping between FE type and element shape
     * @param shape_type
     * @param order
     * @param fe_type
     */
    void addFeType(MeshLib::MeshElemType shape_type, int order, int fe_type)
    {
        MyKey key(shape_type, order);
        _map_shape2fetype[key] = fe_type;
    }

private:
    MyTable _map_shape2fetype;
};

}
