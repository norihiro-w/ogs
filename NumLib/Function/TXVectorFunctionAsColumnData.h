/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 */

#pragma once

#include "TXPosition.h"
#include "ITXFunction.h"

namespace NumLib
{

template <class F_VECTOR>
class TXVectorFunctionAsColumnData : public ITXFunction
{
public:
    explicit TXVectorFunctionAsColumnData(F_VECTOR* f_vec, size_t column_id)
    : _f_vec(f_vec), _column_id(column_id)
    {
        ITXFunction::isConst(f_vec->isConst());
        ITXFunction::isTemporallyConst(f_vec->isTemporallyConst());
        ITXFunction::isSpatiallyConst(f_vec->isSpatiallyConst());
    };

    virtual ~TXVectorFunctionAsColumnData() {};

    void resetVectorFunction(F_VECTOR* f_vec) {_f_vec = f_vec; };

    virtual void eval(const TXPosition x, DataType &v) const
    {
        DataType tmp_v;
        _f_vec->eval(x, tmp_v);
        if (tmp_v.array().size()>0) {
            v.resize(1,1);
            v(0,0) = tmp_v(_column_id);
        }
    }

    virtual TXVectorFunctionAsColumnData<F_VECTOR>* clone() const
    {
        return new TXVectorFunctionAsColumnData<F_VECTOR>(_f_vec, _column_id);
    }

private:
    F_VECTOR* _f_vec;
    size_t _column_id;
};

} //end


