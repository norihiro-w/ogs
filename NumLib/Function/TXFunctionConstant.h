/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file TXFunctionConstant.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include "TXPosition.h"
#include "ITXFunction.h"

namespace NumLib
{

/**
 * \brief Constant value
 */
class TXFunctionConstant : public ITXFunction
{
public:
    /**
     *
     * @param val
     */
    explicit TXFunctionConstant(double val) : _vec(1,1)
    {
        _vec(0,0) = val;
        ITXFunction::isConst(true);
        ITXFunction::isTemporallyConst(true);
        ITXFunction::isSpatiallyConst(true);
    };

    /**
     *
     * @param val
     */
    explicit TXFunctionConstant(const DataType &val) : _vec(val) {};

    /**
     * Copy constructor
     * @param src
     */
    TXFunctionConstant(const TXFunctionConstant &src) : _vec(src._vec) {};

    /**
     *
     */
    virtual ~TXFunctionConstant() {};

    /**
     *
     * @return
     */
    virtual TXFunctionConstant* clone() const
    {
        return new TXFunctionConstant(*this);
    }

    /**
     *
     * @param
     * @param v
     */
    virtual void eval(const TXPosition /*x*/, DataType &v) const
    {
        v = _vec;
    }

    /**
     *
     * @param
     * @param val
     */
    virtual void eval(const TXPosition /*x*/, double &val) const
    {
        val = _vec(0,0);
    }

    /**
     *
     * @param
     * @param val
     */
    virtual void eval(const double* /*x*/, double &val) const
    {
        val = _vec(0,0);
    }

//    /**
//     *
//     * @param v
//     */
//    void eval(double &v) const
//    {
//        v = _vec(0,0);
//    };
//
//    /**
//     *
//     * @param v
//     */
//    void eval(DataType &v) const
//    {
//        v = _vec;
//    };

private:
    DataType _vec;
};

} //end


