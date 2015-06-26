/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file TXFunctionTimeCurve
 *
 * Created on 2012-11-30 by Norihiro Watanabe
 */

#pragma once

#include <string>
#include <vector>
#include "TXPosition.h"
#include "ITXFunction.h"

namespace NumLib
{

/**
 * \brief Time-curve but spatially invariant function
 */
template <class T_INTERPOLATION>
class TXFunctionTimeCurve : public ITXFunction
{
public:
    /**
     *
     * @param interpolate
     */
    explicit TXFunctionTimeCurve(T_INTERPOLATION* interpolate)
    : _interpolate(interpolate)
    {
        ITXFunction::isConst(false);
        ITXFunction::isTemporallyConst(false);
        ITXFunction::isSpatiallyConst(true);
    };

    /**
     *
     * @param src
     */
    TXFunctionTimeCurve(const TXFunctionTimeCurve &src)
    : _interpolate(new T_INTERPOLATION(*src._interpolate))
    {
        ITXFunction::isConst(false);
        ITXFunction::isTemporallyConst(false);
        ITXFunction::isSpatiallyConst(true);
    }

    /**
     *
     */
    virtual ~TXFunctionTimeCurve()
    {
    	delete _interpolate;
    };

    /**
     *
     * @param x
     * @param val
     */
    virtual void eval(const TXPosition x, double &val) const
    {
        val = _interpolate->getValue(x.getTime());
    }

    /**
     *
     * @return
     */
    virtual TXFunctionTimeCurve* clone() const
    {
        return new TXFunctionTimeCurve(*this);
    }

private:
    T_INTERPOLATION* _interpolate;
};

} //end


