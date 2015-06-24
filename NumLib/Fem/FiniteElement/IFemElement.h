/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 */

#pragma once

#include <vector>

#include "MathLib/DataType.h"
#include "NumLib/Fem/CoordinatesMapping/ShapeMatrices.h"

namespace MeshLib
{
class Element;
}

namespace NumLib
{
class IIntegration;

class IFiniteElement
{
public:
	IFiniteElement() : _integration(nullptr) {}
    virtual ~IFiniteElement() {}

    /// return current mesh element
    virtual const MeshLib::Element* getMeshElement() const = 0;

    /// Sets the mesh element
    virtual void setMeshElement(const MeshLib::Element &e) = 0;

    void setIntegrationMethod(IIntegration* integrate) {_integration = integrate;}

    IIntegration* getIntegrationMethod() {return _integration;}

    /**
     * compute shape functions
     *
     * @param natural_pt    position in natural coordinates
     * @param shape         evaluated shape function matrices
     */
    virtual void computeShapeFunctionsd(const double *natural_pt, DynamicShapeMatrices &shape) const = 0;

private:
    IIntegration* _integration;
};

} //end

