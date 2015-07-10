/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef DENSEVECTOR_H_
#define DENSEVECTOR_H_

#include <vector>
#include <valarray>
#include <fstream>
#include <iterator>

#include "MathLib/LinAlg/IVector.h"
#include "MathLib/LinAlg/VectorNorms.h"

namespace MathLib
{

/**
 * Dense vector class
 */
template <typename T>
class DenseVector : public IVector
{
public:
	typedef T FP_T;

public:
	/**
	 * Constructor for initialization of the number of rows
	 * @param nrows number of rows
	 */
	explicit DenseVector(std::size_t nrows=0)
	: _vec(nrows)
	{}

	DenseVector(const DenseVector<T> &src)
	: _vec(src._vec)
	{}

	LinAlgLibType getLinAlgLibType() const override {return LinAlgLibType::Dense;}

	/// duplicate this vector
	IVector* duplicate() const override
	{
		return new DenseVector(*this);
	}

	/// return a start index of the active data range
	std::size_t getRangeBegin() const override { return 0;}

	/// return an end index of the active data range
	std::size_t getRangeEnd() const override { return this->size(); }

	/// get entry
	T get(std::size_t i) const override { return (*this)[i]; }

	/// set a value to entry
	void set(std::size_t i, T v) override { (*this)[i] = v; }

	/// add a value to entry
	void add(std::size_t i, T v) override { (*this)[i] += v; }

	/**
	 * add a sub vector
	 * @param pos       positions of each sub-vector entry in this vector
	 * @param sub_vec   sub-vector
	 */
	template<class T_SUBVEC>
	void add(const std::vector<std::size_t> &pos, const T_SUBVEC &sub_vec)
	{
		for (std::size_t i=0; i<pos.size(); ++i) {
			this->add(pos[i], sub_vec[i]);
		}
	}

	/**
	 * writes the matrix entries into a file
	 * @param filename output file name
	 */
	void write (const std::string &filename) const override
	{
		std::ofstream os(filename);
		//os << *this;
		os.close();
	}

	/// vector operation: add
	void operator+= (const IVector& v) override
	{
		_vec += static_cast<const DenseVector<T>&>(v)._vec;
	}

	/// vector operation: subtract
	void operator-= (const IVector& v) override
	{
		_vec -= static_cast<const DenseVector<T>&>(v)._vec;
	}

	/// set all values in this vector
	IVector& operator= (T v) override { _vec = v; return *this; }

	/// vector operation: set data
	IVector& operator= (const IVector &src) override
	{
		_vec = static_cast<const DenseVector<T>&>(src)._vec;
		return *this;
	}

	/// set all values in this vector
	IVector& operator*= (T v) override { _vec *= v; return *this; }

	/// access entry
	T operator[] (std::size_t rowId) const override {return _vec[rowId];}

	/// access entry
	T& operator[] (std::size_t rowId) {return _vec[rowId];}

	/// returns the dimension of this vector
	std::size_t size() const override {return _vec.size();}

	/// returns a raw data of this vector
	std::valarray<T>& getRawVector() {return _vec;}

	/// returns a raw data of this vector
	const std::valarray<T>& getRawVector() const {return _vec;}

	/// returns norm1
	T norm1() const override { return MathLib::norm_1(&_vec[0], _vec.size()); }

	/// returns norm2
	T norm2() const override { return MathLib::norm_2(&_vec[0], _vec.size()); }

	/// returns maximum norm
	T norm_max() const override { return MathLib::norm_max(&_vec[0], _vec.size()); }

private:
	std::valarray<T> _vec;
};

/// returns norm1
template<class T>
inline double norm_1(const DenseVector<T> &v) { return MathLib::norm_1(&v.getRawVector()[0], v.size()); }

/// returns norm2
template<class T>
inline double norm_2(const DenseVector<T> &v) { return MathLib::norm_2(&v.getRawVector()[0], v.size()); }

/// returns maximum norm
template<class T>
inline double norm_max(const DenseVector<T> &v) { return MathLib::norm_max(&v.getRawVector()[0], v.size()); }

} //MathLib

namespace std
{

template <typename T>
inline double* begin(MathLib::DenseVector<T>& v)
{
	return std::begin(v.getRawVector());
}

template <typename T>
inline double* end(MathLib::DenseVector<T>& v)
{
	return std::end(v.getRawVector());
}

} //std

/**
 * writes a vector content into the output stream
 * @param os the output stream
 */
template <typename T>
std::ostream& operator<<(std::ostream& os, MathLib::DenseVector<T> const & v)
{
	std::copy(std::begin(v), std::end(v), std::ostream_iterator<T>(os, "\n"));
	return os;
}


#endif /* DENSEVECTOR_H_ */
