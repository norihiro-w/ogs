/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef EIGENMATRIX_H_
#define EIGENMATRIX_H_

#include <cassert>
#ifndef NDEBUG
#include <fstream>
#include <string>
#endif

#include <Eigen/Sparse>

#include "MathLib/LinAlg/IMatrix.h"
#include "MathLib/LinAlg/RowColumnIndices.h"
#include "EigenVector.h"

namespace MathLib
{

/**
 * Global matrix based on Eigen sparse matrix
 *
 * The matrix will be dynamically allocated during construction.
 */
class EigenMatrix final : public IMatrix
{
public:
    using RawMatrixType = Eigen::SparseMatrix<double, Eigen::RowMajor>;

    /**
     * constructor
     * @param n the number of rows (that is equal to the number of columns)
     */
    explicit EigenMatrix(std::size_t n) :_mat(n, n) {}

    virtual ~EigenMatrix() {}

    virtual LinAlgLibType getLinAlgLibType() const { return LinAlgLibType::Eigen; }
     
    /// return the number of rows
    std::size_t getNRows() const { return _mat.rows(); }

    /// return the number of columns
    std::size_t getNCols() const { return _mat.cols(); }

    /// return a start index of the active data range
    std::size_t getRangeBegin() const  { return 0; }

    /// return an end index of the active data range
    std::size_t getRangeEnd() const  { return getNRows(); }

    /// reset data entries to zero.
    void setZero()
    {
        _mat.setZero();
    }

    /// set a value to the given entry. If the entry doesn't exist, this class
    /// dynamically allocates it.
    int set(std::size_t row, std::size_t col, double val)
    {
        _mat.coeffRef(row, col) = val;
        return 0;
    }

    /// add a value to the given entry. If the entry doesn't exist, the value is
    /// inserted.
    int add(std::size_t row, std::size_t col, double val)
    {
        assert(row < getNRows() && col < getNCols());
        _mat.coeffRef(row, col) += val;
        return 0;
    }

    /// get value. This function returns zero if the element doesn't exist.
    double get(std::size_t row, std::size_t col) const
    {
        assert(row < getNRows() && col < getNCols());
        return _mat.coeff(row, col);
    }

    /// get value. This function returns zero if the element doesn't exist.
    double operator() (std::size_t row, std::size_t col) const
    {
        return get(row, col);
    }

    /// get a maximum value in diagonal entries
    double getMaxDiagCoeff()
    {
        return _mat.diagonal().maxCoeff();
    }

    /// y = mat * x
    void multiply(const IVector &x, IVector &y) const
    {
        auto &xx = static_cast<const EigenVector&>(x);
        auto &yy = static_cast<EigenVector&>(y);
        yy.getRawVector() = _mat * xx.getRawVector();
    }

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

    void zeroRows(const std::vector<std::size_t> &vec_knownX_id)
    {
        for (auto rowId : vec_knownX_id)
        {
            for (Eigen::SparseMatrix<double>::InnerIterator it(_mat,rowId); it; ++it)
                if (it.col()!=rowId)
                    it.valueRef() = 0.0;
        }
    }

    /// return always true, i.e. the matrix is always ready for use
    bool isAssembled() const { return true; }

    /// assemble the matrix
    virtual void assemble() {}

//#ifndef NDEBUG
    /// printout this matrix for debugging
    void write(const std::string &filename) const
    {
        std::ofstream of(filename);
        if (of)
            write(of);
    }

    /// printout this matrix for debugging
    void write(std::ostream &os) const
    {
        for (int k=0; k<_mat.outerSize(); ++k)
          for (Eigen::SparseMatrix<double>::InnerIterator it(_mat,k); it; ++it)
              os << it.row() << " " << it.col() << ": " << it.value() << "\n";
        os << std::endl;
    }
//#endif

    RawMatrixType& getRawMatrix() { return _mat; }
    const RawMatrixType& getRawMatrix() const { return _mat; }

protected:
    RawMatrixType _mat;
};

template<class T_DENSE_MATRIX>
void
EigenMatrix::add(std::vector<std::size_t> const& row_pos, std::vector<std::size_t> const& col_pos,
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


} // end namespace MathLib

#endif

