/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file TXFunction.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include "MathLib/DataType.h"
#include "NumLib/Function/IFunction.h"
#include "TXPosition.h"

namespace NumLib
{

/**
 * \brief Interface of any functions in space-time domain
 *
 * This class aims to be an abstract of spatially and temporally distributed data such as
 * physical quantity (e.g. head, stress) and material property (e.g. permeability).
 *
 * TXFunction
 * - is evaluated at particular position in space-time domain
 * - returns scalar or vector value
 * - has some attributes (e.g constant)
 *
 */
class ITXFunction : public IClonable
{
public:
    typedef MathLib::LocalMatrix DataType;
    
    ///
    ITXFunction() : _is_const(false), _is_temporally_const(false), _is_spatially_const(false) {}

    ///
    virtual ~ITXFunction() {}

    /// evaluate this function at the given position and return vector data
    ///
    /// \param x  position in space and time
    /// \param v  evaluated vector
    virtual void eval(const TXPosition /*x*/, DataType &/*v*/) const {}

    /// evaluate this function at the given position and return scalar value.
    ///
    /// Default behavior of this function is to return the 1st component of the vector.
    /// \param x  position in space and time
    /// \param v  evaluated scalar
    virtual void eval(const TXPosition x, double &v) const
    {
        DataType tmp(1,1);
        this->eval(x, tmp);
#ifdef OGS_USE_EIGEN
        v = tmp(0);
#else
        v = tmp(0,0);
#endif
    }

    ///
    virtual void eval(const double* x, DataType &v) const
    {
        TXPosition pos(x);
        this->eval(pos, v);
    }

    ///
    virtual void eval(const double* x, double &v) const
    {
        TXPosition pos(x);
        this->eval(pos, v);
    }

    ///
    virtual ITXFunction* clone() const = 0;

    ///
    bool isConst() const {return _is_const;}
    void isConst(bool b) {_is_const = b;}
    ///
    bool isTemporallyConst() const {return _is_temporally_const;}
    void isTemporallyConst(bool b) {_is_temporally_const = b;}
    ///
    bool isSpatiallyConst() const {return _is_spatially_const;}
    void isSpatiallyConst(bool b) {_is_spatially_const = b;}

private:
    bool _is_const;
    bool _is_temporally_const;
    bool _is_spatially_const;
};

} //end


