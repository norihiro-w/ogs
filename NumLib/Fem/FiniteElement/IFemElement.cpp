/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 */

#include "IFemElement.h"

#include "MeshLib/Elements/Element.h"
#include "NumLib/Fem/Integration/GaussIntegrationPolicy.h"

namespace NumLib
{

IFiniteElement::~IFiniteElement()
{
    delete _integration;
}

IIntegration& IFiniteElement::getIntegrationMethod() const
{
    if (_integration == nullptr)
        _integration = NumLib::getIntegrationMethod(getMeshElement()->getCellType());
    return *_integration;
}

} //end

