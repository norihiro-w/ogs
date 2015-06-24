/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 */

#pragma once

namespace MeshLib
{
class Element;
}

namespace NumLib
{

class IMeshElementToFemElementType
{
public:
    virtual ~IMeshElementToFemElementType() {};

    virtual int getFeType(const MeshLib::Element& e, const int order) const = 0;
};


}
