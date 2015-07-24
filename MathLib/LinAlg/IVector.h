/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef IVECTOR_H_
#define IVECTOR_H_

#include <iostream>
#include <vector>

#include "BaseLib/CodingTools.h"
#include "MathLib/DataType.h"
#include "LinAlgLibType.h"

namespace MathLib
{

class IVector
{
public:
    virtual ~IVector() {}

    virtual LinAlgLibType getLinAlgLibType() const = 0;

    virtual IVector* duplicate() const = 0;

    /// return a vector length
    virtual std::size_t size() const = 0;

    /// return a start index of the active data range
    virtual std::size_t getRangeBegin() const = 0;

    /// return an end index of the active data range
    virtual std::size_t getRangeEnd() const = 0;

    /// set all values in this vector
    virtual IVector& operator= (double v) = 0;

    /// set all values in this vector
    virtual IVector& operator*= (double v) = 0;

    /// access entry
    virtual double operator[] (std::size_t rowId) const = 0;

    /// get entry
    virtual double get(std::size_t rowId) const = 0;

    /// set entry
    virtual void set(std::size_t rowId, double v) = 0;

    /// add entry
    virtual  void add(std::size_t rowId, double v) = 0;

    /// printout this equation for debugging
    virtual void write (const std::string &filename) const = 0;

    /// vector operation: set data
    virtual IVector& operator= (const IVector &/*src*/) {return*this;} //TODO clang bug? no need to have definitio

    /// vector operation: add
    virtual void operator+= (const IVector& v) = 0;

    /// vector operation: subtract
    virtual void operator-= (const IVector& v) = 0;

    ///
    template<class T_SUBVEC>
    void add(const std::vector<std::size_t> &pos, const T_SUBVEC &sub_vec, double fac = 1.0)
    {
        for (std::size_t i=0; i<pos.size(); ++i) {
            this->add(pos[i], fac * sub_vec[i]);
        }
    }

    virtual double norm1() const = 0;
    virtual double norm2() const = 0;
    virtual double norm_max() const = 0;

    virtual void assemble() = 0;
};

} // MathLib

#endif //IVector

