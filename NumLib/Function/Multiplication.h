/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file Multiplication.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

namespace NumLib
{

/**
 *
 */
template <typename Tval>
class Multiplication
{
public:
    void doit(const Tval &v1, const Tval &v2, Tval &val)
    {
        val = v1 * v2;
    }
};

} //end


