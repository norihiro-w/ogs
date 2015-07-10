/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef IMATRIX_H_
#define IMATRIX_H_

#include <iostream>
#include <cmath>
#include <vector>

#include "RowColumnIndices.h"
#include "IVector.h"

namespace MathLib
{

class IMatrix
{
public:
    virtual ~IMatrix(){}

    virtual LinAlgLibType getLinAlgLibType() const = 0;

    /// return the number of rows
    virtual std::size_t getNRows() const = 0;

    /// return the number of columns
    virtual std::size_t getNCols() const = 0;

    /// return a start index of the active data range
    virtual std::size_t getRangeBegin() const = 0;

    /// return an end index of the active data range
    virtual std::size_t getRangeEnd() const = 0;

    /// reset this matrix with keeping its original dimension
    virtual void setZero() = 0;

    /// set entry
    virtual int set(std::size_t rowId, std::size_t colId, double v) = 0;

    /// add value
    virtual int add(std::size_t rowId, std::size_t colId, double v) = 0;

    /// printout this equation for debugging
    virtual void write(const std::string &filename) const = 0;

    /// get a maximum value in diagonal entries
    virtual double getMaxDiagCoeff() = 0;

    /// y = mat * x
    virtual void multiply(const IVector &x, IVector &y) const = 0;

    /// Add sub-matrix at positions \c row_pos and same column positions as the
    /// given row positions.
    template<class T_DENSE_MATRIX>
    void add(std::vector<std::size_t> const& row_pos,
            const T_DENSE_MATRIX &sub_matrix,
            double fkt = 1.0)
    {
        this->add(row_pos, row_pos, sub_matrix, fkt);
    }

    /// Add sub-matrix at positions given by \c indices.
    template<class T_DENSE_MATRIX>
    void add(RowColumnIndices<std::size_t> const& indices,
            const T_DENSE_MATRIX &sub_matrix,
            double fkt = 1.0)
    {
        this->add(indices.rows, indices.columns, sub_matrix, fkt);
    }

    ///
    template <class T_DENSE_MATRIX>
    void add(std::vector<std::size_t> const& row_pos,
            std::vector<std::size_t> const& col_pos, const T_DENSE_MATRIX &sub_matrix,
            double fkt = 1.0);

    /// return if this matrix is already assembled or not
    virtual bool isAssembled() const = 0;

    /// assemble the matrix
    virtual void assemble() = 0;

    /// zeros all entries of a set of rows of a matrix
    virtual void zeroRows(const std::vector<std::size_t> &vec_knownX_id) = 0;
};

template<class T_DENSE_MATRIX>
void
IMatrix::add(std::vector<std::size_t> const& row_pos, std::vector<std::size_t> const& col_pos,
        const T_DENSE_MATRIX &sub_matrix, double fkt)
{
    const std::size_t n_rows = row_pos.size();
    const std::size_t n_cols = col_pos.size();
    for (std::size_t i = 0; i < n_rows; i++) {
        const std::size_t row = row_pos[i];
        for (std::size_t j = 0; j < n_cols; j++) {
            const std::size_t col = col_pos[j];
            add(row, col, fkt * sub_matrix(i, j));
        }
    }
};


} // MathLib

#endif //IMATRIX_H_

