/**
 * \file
 * \author Norihiro Watanabe
 * \date   2013-04-16
 * \brief
 *
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <gtest/gtest.h>

#ifndef TESTTOOLS_H_
#define TESTTOOLS_H_

#define ASSERT_ARRAY_NEAR(E,A,N,eps)\
    for (size_t i=0; i<(unsigned)(N); i++) \
        ASSERT_NEAR((E)[i], (A)[i], (eps));

#define ASSERT_ARRAY_EQ(E,A,N)\
    for (size_t i=0; i<(unsigned)(N); i++) \
        ASSERT_EQ((E)[i], (A)[i]);

template <class T>
std::vector<double> to_array(const T &m)
{
#ifdef OGS_USE_EIGEN
    auto cols = m.cols();
#else
    auto cols = m.columns();
#endif
    std::vector<double> v(m.rows()*cols);
    for (unsigned i=0;i<m.rows();i++)
        for (unsigned j=0;j<cols;j++)
            v[i*cols+j] = m(i,j);
    return v;
}

#endif // TESTTOOLS_H_
