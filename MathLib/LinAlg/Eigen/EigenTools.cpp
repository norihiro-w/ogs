/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "EigenTools.h"

#include <logog/include/logog.hpp>

#include "EigenMatrix.h"
#include "EigenVector.h"

namespace MathLib
{

void applyKnownSolution(EigenMatrix &A_, EigenVector &b_, const std::vector<std::size_t> &vec_knownX_id,
		const std::vector<double> &vec_knownX_x, double /*penalty_scaling*/)
{
    using SpMat = EigenMatrix::RawMatrixType;
    auto &A = A_.getRawMatrix();
    auto &b = b_.getRawVector();
    const std::size_t n_rows = A.rows();
    for (std::size_t ix=0; ix<vec_knownX_id.size(); ix++)
    {
        auto row_id = vec_knownX_id[ix];
        auto x = vec_knownX_x[ix];
        //A(k, j) = 0.
        for (SpMat::InnerIterator it(A,row_id); it; ++it)
            it.valueRef() = .0;
        //b_i -= A(i,k)*val, i!=k
        for (std::size_t i=0; i<n_rows; i++)
            for (SpMat::InnerIterator it(A,i); it; ++it)
            {
                if (it.col()!=row_id) continue;
                b[i] -= it.value()*x;
            }
        //b_k = val
        b[row_id] = x;
        //A(i, k) = 0., i!=k
        for (std::size_t i=0; i<n_rows; i++)
            for (SpMat::InnerIterator it(A,i); it; ++it)
            {
                if (it.col()!=row_id) continue;
                it.valueRef() = 0.0;
            }
        //A(k, k) = 1.0
        A.coeffRef(row_id, row_id) = 1.0; //=x
    }
}

void buildCRSMatrixFromEigenMatrix(const EigenMatrix &A_, int &nonzero, int*& row_ptr, int*& col_idx, double*& data)
{
    auto &A = A_.getRawMatrix();
    const size_t dimension = A.rows();

    for (size_t i=0; i<A.outerSize(); ++i)
      for (EigenMatrix::RawMatrixType::InnerIterator it(A,i); it; ++it)
        nonzero++;

    data = new double [nonzero];
    size_t temp_cnt = 0;
    for (size_t i=0; i<A.outerSize(); ++i)
      for (EigenMatrix::RawMatrixType::InnerIterator it(A,i); it; ++it)
        data[temp_cnt++] = it.value();

    assert(temp_cnt==nonzero);
    row_ptr = new int[A.rows()+1];
    col_idx = new int[nonzero];
    long counter_ptr = 0, counter_col_idx = 0;
    long cnt_row = 0;

    for (size_t i=0; i<A.outerSize(); ++i)
    {
        row_ptr[cnt_row++] = counter_ptr;         // starting point of the row
      for (EigenMatrix::RawMatrixType::InnerIterator it(A,i); it; ++it)
      {
        col_idx[counter_col_idx] = it.col();
        ++counter_ptr;
        ++counter_col_idx;
      }
    }
    row_ptr[A.rows()] = counter_ptr;
#if 0
    //output CRS
    cout << "PTR:" << endl;
    for (size_t i=0; i<A.rows()+1; i++)
        cout << ptr[i] << ", ";
    cout << endl;
    cout << "ColID:" << endl;
    for (size_t i=0; i<nonzero; i++)
        cout << col_idx[i] << ", ";
    cout << endl;
    cout << "Data:" << endl;
    for (size_t i=0; i<nonzero; i++)
        cout << crs_data[i] << ", ";
    cout << endl;
#endif
}

} // MathLib



