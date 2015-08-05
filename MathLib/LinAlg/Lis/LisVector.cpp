/**
 * \file
 * \author Norihiro Watanabe
 * \date   2013-04-16
 * \brief  Implementation of the LisVector class.
 *
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "LisVector.h"
#include "LisCheck.h"

namespace MathLib
{

LisVector::LisVector(std::size_t length)
{
	lis_vector_create(0, &_vec);
	lis_vector_set_size(_vec, 0, length);
}

LisVector::LisVector(std::size_t length, double* data)
{
    lis_vector_create(0, &_vec);
    lis_vector_set_size(_vec, 0, length);
    for (std::size_t i=0; i<length; i++)
        lis_vector_set_value(LIS_INS_VALUE, i, data[i], _vec);
}

LisVector::LisVector(LisVector const &src)
{
	lis_vector_duplicate(src._vec, &_vec);
	lis_vector_copy(src._vec, _vec);
}

LisVector::~LisVector()
{
	lis_vector_destroy(_vec);
}

IVector* LisVector::duplicate() const
{
	return new LisVector(*this);
}

IVector& LisVector::operator= (const IVector &src)
{
	lis_vector_copy(static_cast<const LisVector&>(src)._vec, _vec);
	return *this;
}

void LisVector::operator+= (const IVector& v)
{
	lis_vector_axpy(1.0, static_cast<const LisVector&>(v)._vec, _vec);
}

void LisVector::operator-= (const IVector& v)
{
	lis_vector_axpy(-1.0, static_cast<const LisVector&>(v)._vec, _vec);
}

IVector& LisVector::operator= (double v)
{
	lis_vector_set_all(v, _vec);
	return *this;
}

IVector& LisVector::operator*= (double v)
{
	lis_vector_scale(v, _vec);
	return *this;
}

std::size_t LisVector::size() const
{
	LIS_INT dummy;
	LIS_INT size;
	int const ierr = lis_vector_get_size(_vec, &dummy, &size);
	checkLisError(ierr);
	return size;
}

void LisVector::write (const std::string &filename) const
{
	lis_output_vector(_vec, LIS_FMT_PLAIN, const_cast<char*>(filename.c_str()));
}


} // MathLib

