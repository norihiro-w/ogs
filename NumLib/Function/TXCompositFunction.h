/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file TXCompositFunction.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include "TXPosition.h"
#include "ITXFunction.h"

namespace NumLib
{

/**
 * \brief Temporal-spatial function made of two functions
 *
 */
template <class T1, class T2, template <typename> class T_OPERATOR>
class TXCompositFunction : public ITXFunction
{
public:
    /**
     *
     * @param f1
     * @param f2
     */
    TXCompositFunction(T1* f1, T2* f2) : _f1(f1->clone()), _f2(f2->clone())
    {
        isConst(_f1->isConst() && _f2->isConst());
        isTemporallyConst(_f1->isTemporallyConst() && _f2->isTemporallyConst());
        isSpatiallyConst(_f1->isSpatiallyConst() && _f2->isSpatiallyConst());
    };

    /**
     * Copy constructor
     * @param src
     */
    TXCompositFunction(const TXCompositFunction &src)
    : _f1(src._f1->clone()), _f2(src._f2->clone())
    {
        isConst(_f1->isConst() && _f2->isConst());
        isTemporallyConst(_f1->isTemporallyConst() && _f2->isTemporallyConst());
        isSpatiallyConst(_f1->isSpatiallyConst() && _f2->isSpatiallyConst());
    }

    /**
     *
     */
    virtual ~TXCompositFunction()
    {
    	delete _f1;
    	delete _f2;
    };

    /**
     *
     * @param x
     * @param val
     */
    virtual void eval(const TXPosition x, DataType &val) const
    {
        DataType v1, v2;
        _f1->eval(x, v1);
        _f2->eval(x, v2);
        T_OPERATOR<DataType> op;
        op.doit(v1, v2, val);
    }

    /**
     *
     * @param x
     * @param v
     */
    virtual void eval(const TXPosition x, double &val) const
    {
        double v1, v2;
        _f1->eval(x, v1);
        _f2->eval(x, v2);
        T_OPERATOR<double> op;
        op.doit(v1, v2, val);
    }

    /**
     *
     * @param x
     * @param val
     */
    virtual void eval(const double* x, double &val) const
    {
        double v1, v2;
        _f1->eval(x, v1);
        _f2->eval(x, v2);
        T_OPERATOR<double> op;
        op.doit(v1, v2, val);
    }

    /**
     *
     * @return
     */
    virtual TXCompositFunction* clone() const
    {
        return new TXCompositFunction<T1, T2, T_OPERATOR>(*this);
    }

private:
    T1* _f1;
    T2* _f2;
};

} //end


