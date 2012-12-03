/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file TXFunctionDirect.h
 *
 * Created on 2012-09-22 by Norihiro Watanabe
 */
 
 #pragma once

#include "DiscreteLib/Core/IDiscreteVector.h"
#include "TXPosition.h"
#include "ITXDiscreteFunction.h"

namespace NumLib
{

template <class T>
class TXFunctionDirect : public ITXDiscreteFunction<T>
{
public:
    typedef DiscreteLib::IDiscreteVector<T> MyVector;
    typedef MathLib::LocalMatrix DataType;
    
    explicit TXFunctionDirect(const MyVector* direct_data)
    : _direct_data(direct_data)
    {
        ITXFunction::isConst(true);
        ITXFunction::isTemporallyConst(true);
        ITXFunction::isSpatiallyConst(false);
    };

    virtual ~TXFunctionDirect() {};

    virtual void eval(const TXPosition x, DataType &v) const
    {
        switch (x.getIdObjectType()) {
        case NumLib::TXPosition::Node:
        case NumLib::TXPosition::Element:
            {
                v = (*_direct_data)[x.getId()];
            }
            break;
        default:
            break;
        }
    }
    
    virtual TXFunctionDirect<T>* clone() const
    {
        return new TXFunctionDirect<T>(_direct_data);
    }

    virtual MyVector* getDiscreteData() 
    {
        return (MyVector*)_direct_data;
    }
    
    virtual const MyVector* getDiscreteData() const
    {
        return _direct_data;
    }
    
private:
    const MyVector* _direct_data;
};

template <>
class TXFunctionDirect<double> : public ITXDiscreteFunction<double>
{
public:
    typedef DiscreteLib::IDiscreteVector<double> MyVector;
    
    explicit TXFunctionDirect(const MyVector* direct_data)
    : _direct_data(direct_data)
    {
        ITXFunction::isConst(true);
        ITXFunction::isTemporallyConst(true);
        ITXFunction::isSpatiallyConst(false);
    };

    virtual ~TXFunctionDirect() {};

    virtual void eval(const TXPosition x, DataType &v) const
    {
        switch (x.getIdObjectType()) {
        case NumLib::TXPosition::Node:
        case NumLib::TXPosition::Element:
            {
                v(0,0) = (*_direct_data)[x.getId()];
            }
            break;
        default:
            break;
        }
    }
    
    virtual TXFunctionDirect<double>* clone() const
    {
        return new TXFunctionDirect<double>(_direct_data);
    }

    virtual MyVector* getDiscreteData() 
    {
        return (MyVector*)_direct_data;
    }
    
    virtual const MyVector* getDiscreteData() const
    {
        return _direct_data;
    }
    
private:
    const MyVector* _direct_data;
};

} //end
