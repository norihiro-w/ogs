/**
 * \file
 * \author Norihiro Watanabe
 * \date   2013-04-16
 * \brief  Definition of the LisVector class.
 *
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef LISVECTOR_H_
#define LISVECTOR_H_

#include <string>
#include <vector>

#include <lis.h>

#include "MathLib/LinAlg/IVector.h"

namespace MathLib
{

/**
 * \brief Lis vector wrapper class
 */
class LisVector : public IVector
{
public:
	/**
	 * Constructor for initialization of the number of rows
	 * @param length number of rows
	 */
    explicit LisVector(std::size_t length);

    /**
     * Constructor using the given raw data
     * @param length the length of the vector
     * @param data   the raw data
     */
    LisVector(std::size_t length, double* data);

    /// copy constructor
    LisVector(LisVector const &src);

    /**
     *
     */
    virtual ~LisVector();

    virtual LinAlgLibType getLinAlgLibType() const {return LinAlgLibType::Lis;}

    /// duplicate this vector
    IVector* duplicate() const;

    /// return a vector length
    std::size_t size() const;

    /// return a start index of the active data range
    std::size_t getRangeBegin() const { return 0;}

    /// return an end index of the active data range
    std::size_t getRangeEnd() const { return this->size(); }

    /// set all values in this vector
    IVector& operator= (double v);

    /// set all values in this vector
    IVector& operator*= (double v);

    /// access entry
    double operator[] (std::size_t rowId) const { return get(rowId); }

    /// get entry
    double get(std::size_t rowId) const
    {
        double v = .0;
        lis_vector_get_value(_vec, rowId, &v);
        return v;
    }

    /// set entry
    void set(std::size_t rowId, double v)
    {
        lis_vector_set_value(LIS_INS_VALUE, rowId, v, _vec);
    }

    /// add entry
    void add(std::size_t rowId, double v)
    {
        lis_vector_set_value(LIS_ADD_VALUE, rowId, v, _vec);
    }

    /// printout this equation for debugging
    void write (const std::string &filename) const;

    /// return a raw Lis vector object
    LIS_VECTOR& getRawVector() {return _vec; }

    /// vector operation: set data
    IVector& operator= (const IVector &src);

    /// vector operation: add
    void operator+= (const IVector& v);

    /// vector operation: subtract
    void operator-= (const IVector& v);

    ///
    template<class T_SUBVEC>
    void add(const std::vector<std::size_t> &pos, const T_SUBVEC &sub_vec)
    {
        for (std::size_t i=0; i<pos.size(); ++i) {
            this->add(pos[i], sub_vec[i]);
        }
    }

    void assemble() {}

    double norm1() const
    {
        double n = .0;
        lis_vector_nrm1(_vec, &n);
        return n;
    }
    double norm2() const
    {
        double n = .0;
        lis_vector_nrm2(_vec, &n);
        return n;
    }
    double norm_max() const
    {
        double n = .0;
        lis_vector_nrmi(_vec, &n);
        return n;
    }

private:
    LIS_VECTOR _vec;
};


inline double norm_1(const LisVector &v)
{
	return v.norm1();
}

inline double norm_2(const LisVector &v)
{
	return v.norm2();
}

inline double norm_max(const LisVector &v)
{
	return v.norm_max();
}

} // MathLib

#endif //LISVECTOR_H_

