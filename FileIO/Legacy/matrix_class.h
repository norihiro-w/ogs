/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file matrix_class.h
 *
 * Created on 2004-08-22 by Wenqing Wang
 */

/*========================================================================
   GeoSys - class Matrix, Sparse matrix (Declaration)
   Task:       Matrix object for full matrices.  If the size of matrix is
             small, this class can do efficient matrix operation.
   Function:   See the declaration below
   Design and programm WW
   03/2010 some improvements TF
   ==========================================================================*/
#ifndef matrix_class_INC
#define matrix_class_INC

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <vector>

namespace ogs5
{

namespace Math_Group
{

class Matrix
{
public:
    Matrix(size_t rows, size_t cols = 1);
    Matrix();
    explicit Matrix(const Matrix& m);
    //
    void resize(size_t rows, size_t cols = 1);
    //
    virtual ~Matrix();
    void ReleaseMemory();                 //06.2010. WW

    // Operators
    virtual void operator= (double a);
    virtual void operator*= (double a);
    virtual void operator/= (double a);
    virtual void operator+= (double a);
    void operator= (const Matrix& m);
    void operator+= (const Matrix& m);
    void operator-= (const Matrix& m);

    void GetTranspose(Matrix& m);

    // vec_result = This*vec. vec_result must be initialized
    void multi(const double* vec, double* vec_result, double fac = 1.0);
    // m_result = this*m. m_result must be initialized
    void multi(const Matrix& m, Matrix& m_result, double fac = 1.0);
    // m_result = this*m1*m2. m_result must be initialized
    void multi(const Matrix& m1, const Matrix& m2, Matrix& m_result);

    // Access to members
    virtual double& operator() (size_t i, size_t j = 0) const;
    void LimitSize(size_t nRows, size_t nCols = 1);

    size_t Rows() const {return nrows; }
    size_t Cols() const {return ncols; }
    size_t Size() const {return size; }

    // Print
    void Write(std::ostream& os = std::cout);
    void Write_BIN(std::fstream& os);
    void Read_BIN(std::fstream& is);
protected:
    size_t nrows, nrows0;
    size_t ncols, ncols0;
    size_t size;
    double* data;
    bool Sym;
};

// Symmetrical matrix. 12-01-2005. WW
class SymMatrix : public Matrix
{
public:
    SymMatrix(size_t dim);
    SymMatrix();
    explicit SymMatrix(const SymMatrix& m);

    void resize(size_t dim);
    ~SymMatrix() {}

    // Operators
    void operator= (double a);
    void operator*= (double a);
    void operator+= (double a);
    void operator= (const SymMatrix& m);
    void operator+= (const SymMatrix& m);
    void operator-= (const SymMatrix& m);
    void LimitSize(size_t dim);

    // Access to members
    double& operator() (size_t i, size_t j) const;
};

class DiagonalMatrix : public Matrix
{
private:
    mutable double dummy_zero;
public:
    DiagonalMatrix(size_t dim);
    DiagonalMatrix();
    explicit DiagonalMatrix(const DiagonalMatrix& m);

    void resize(size_t dim);

    ~DiagonalMatrix() {}

    // Operators
    void operator = (double a);
    void operator *= (double a);
    void operator += (double a);
    void operator = (const DiagonalMatrix& m);
    void operator += (const DiagonalMatrix& m);
    void operator -= (const DiagonalMatrix& m);
    void LimitSize(size_t dim);

    // Access to members
    double& operator() (size_t i, size_t j) const;
    double& operator() (size_t i) const;
};

typedef Matrix Vec;


// End of class Matrix
}
}

//==========================================================================
#endif
