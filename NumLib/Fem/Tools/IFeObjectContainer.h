/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file IFeObjectContainer.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include "NumLib/Fem/FiniteElement/IFemElement.h"

namespace NumLib
{

/**
 * \brief Interface of finite element object containers
 */
class IFeObjectContainer
{
public:
    virtual ~IFeObjectContainer() {};
    /// get a finite element object for the given mesh element
    virtual IFiniteElement* getFeObject(const MeshLib::Element &e) const = 0;
    /// return a clone of this object
    virtual IFeObjectContainer* clone() const = 0;
};

}
