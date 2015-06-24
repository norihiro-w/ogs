/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef IINTEGRATION_H_
#define IINTEGRATION_H_

#include <cmath>
#include <array>

#include "MathLib/TemplateWeightedPoint.h"


namespace NumLib
{

class IIntegration
{
    typedef typename MathLib::TemplateWeightedPoint<double, double, 3>
        WeightedPoint;
public:
    virtual ~IIntegration(){}

    /// Change the integration order.
    virtual void setIntegrationOrder(std::size_t order) = 0;

    /// return current integration order.
    virtual std::size_t getIntegrationOrder() const = 0;

    /// return the number of sampling points
    virtual std::size_t getNPoints() const = 0;

    /// Get coordinates of the integration point.
    ///
    /// @param igp       The integration point index
    /// @return a weighted point
    virtual WeightedPoint
    getWeightedPoint(std::size_t igp) const = 0;

};

} // NumLib


#endif //IINTEGRATION_H_
