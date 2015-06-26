/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 */

#pragma once

namespace NumLib
{

class IClonable
{
public:
    virtual ~IClonable() {}
    virtual IClonable* clone() const = 0;
};

}
