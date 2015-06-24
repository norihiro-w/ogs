/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file matrix_class.cpp
 *
 * Created on 2004-08-22 by Wenqing Wang
 */

/*========================================================================
   GeoSys - class Matrix (Definition)
          class vec
   Task:       Matrix object for full matrices.  If the size of matrix is
             small, this class can do efficient matrix operation.
   Function:   See the definition below
   programming:
   22/08/2004  WW
   ==========================================================================*/

/// Matrix
#include <cfloat>
#include <cmath>
#include <iomanip>
//
//#include "mathlib.h"
#include "matrix_class.h"
//

namespace ogs5
{
namespace Math_Group
{
// Constructors
Matrix::Matrix(size_t rows, size_t cols) :
    nrows (rows), nrows0 (rows), ncols (cols), ncols0 (cols),
    size (nrows * ncols), data (new double[size]), Sym (false)
{
    for(size_t i = 0; i < size; i++)
        data[i] = 0.0;
}

Matrix::Matrix() :
    nrows (0), nrows0 (0), ncols (0), ncols0 (0),
    size (nrows * ncols), data (NULL), Sym (false)
{}

Matrix::Matrix(const Matrix& m) :
    nrows (m.nrows), nrows0 (m.nrows), ncols (m.ncols), ncols0 (m.ncols),
    size (nrows * ncols), data (new double[size]), Sym (m.Sym)
{
    for(size_t i = 0; i < size; i++)
        data[i] = m.data[i];
}

void Matrix::resize(size_t rows, size_t cols)
{
    if (size > 0)
    {
        delete[] data;
        data = NULL;
    }

    Sym = false;
    nrows = rows;
    ncols = cols;
    nrows0 = rows;
    ncols0 = ncols;
    size = nrows * ncols;
    data = new double[size];
    for(size_t i = 0; i < size; i++)
        data[i] = 0.0;
}

Matrix::~Matrix()
{
    delete [] data;
    data = NULL;
}
// 06.2010. WW
void Matrix::ReleaseMemory()
{
    delete [] data;
    data = NULL;
}

void Matrix::operator= (double a)
{
    for(size_t i = 0; i < size; i++)
        data[i] = a;
}
void Matrix::operator *= (double a)
{
    for(size_t i = 0; i < size; i++)
        data[i] *= a;
}
void Matrix::operator /= (double a)
{
    for(size_t i = 0; i < size; i++)
        data[i] /= a;
}
void Matrix::operator += (double a)
{
    for(size_t i = 0; i < size; i++)
        data[i] += a;
}
//
void Matrix::operator = (const Matrix& m)
{
#ifdef gDEBUG
    if(nrows != m.Rows() || ncols != m.Cols())
    {
        std::cout << "\n The sizes of the two matrices are not matched" << std::endl;
        abort();
    }
#endif
    for(size_t i = 0; i < nrows; i++)
        for(size_t j = 0; j < ncols; j++)
            data[i * ncols + j] = m(i,j);
}

//
void Matrix::operator += (const Matrix& m)
{
#ifdef gDEBUG
    if(nrows != m.Rows() || ncols != m.Cols())
    {
        std::cout << "\n The sizes of the two matrices are not matched" << std::endl;
        abort();
    }
#endif
    for(size_t i = 0; i < nrows; i++)
        for(size_t j = 0; j < ncols; j++)
            data[i * ncols + j] += m(i,j);
}

//
void Matrix::operator -= (const Matrix& m)
{
#ifdef gDEBUG
    if(nrows != m.Rows() || ncols != m.Cols()) //Assertion, will be removed
    {
        std::cout << "\n The sizes of the two matrices are not matched" << std::endl;
        abort();
    }
#endif
    for(size_t i = 0; i < nrows; i++)
        for(size_t j = 0; j < ncols; j++)
            data[i * ncols + j] -= m(i,j);
}
//
void Matrix::GetTranspose(Matrix& m)
{
#ifdef gDEBUG
    if(ncols != m.Rows() && nrows != m.Cols())
    {
        std::cout << "\n The sizes of the two matrices are not matched" << std::endl;
        abort();
    }
#endif

    for(size_t i = 0; i < m.Rows(); i++)
        for(size_t j = 0; j < m.Cols(); j++)
            //          m(i,j) = data[j*ncols+i];
            m(i,j) = (*this)(j,i);
}
//
// m_results = this*m. m_results must be initialized
void Matrix::multi(const Matrix& m, Matrix& m_result, double fac)
{
#ifdef gDEBUG
    if(ncols != m.Rows() && nrows != m_result.Rows() && m.Cols() != m_result.Cols())
    {
        std::cout << "\n The sizes of the two matrices are not matched" << std::endl;
        abort();
    }
#endif
    for(size_t i = 0; i < m_result.Rows(); i++)
        for(size_t j = 0; j < m_result.Cols(); j++)
        {
            if(Sym && (j > i))
                continue;
            // m_result(i,j) = 0.0;
            for(size_t k = 0; k < ncols; k++)
                //            m_result(i,j) += fac*data[i*ncols+k]*m(k,j);
                m_result(i,j) += fac * (*this)(i,k) * m(k,j);
        }
}

//
// m_results = this*m1*m2. m_results must be  initialized
void Matrix::multi(const Matrix& m1, const Matrix& m2, Matrix& m_result)
{
#ifdef gDEBUG
    if(ncols != m1.Rows() && m1.Cols() != m2.Rows()
       && m2.Cols() != m_result.Cols() && nrows != m_result.Rows())
    {
        std::cout << "\n The sizes of the two matrices are not matched" << std::endl;
        abort();
    }
#endif
    for(size_t i = 0; i < m_result.Rows(); i++)
        for(size_t j = 0; j < m_result.Cols(); j++)
        {
            if(Sym && (j > i))
                continue;
            //m_result(i,j) = 0.0;
            for(size_t k = 0; k < ncols; k++)
                for(size_t l = 0; l < m2.Rows(); l++)
                    //                m_result(i,j) += data[i*ncols+k]*m1(k,l)*m2(l,j);
                    m_result(i,j) += (*this)(i,k) * m1(k,l) * m2(l,j);
        }
}
// vec_result = This*vec. vec_result must be  initialized
void Matrix::multi(const double* vec, double* vec_result, double fac)
{
    for(int i = 0; (size_t)i < nrows; i++)
        for(int j = 0; (size_t)j < ncols; j++)
            vec_result[i] += fac * (*this)(i,j) * vec[j];
}

double& Matrix::operator() (size_t i, size_t j) const
{
#ifdef gDEBUG
    if(i >= nrows || j >= ncols)
    {
        std::cout << "\n Index exceeds the size of the matrix" << std::endl;
        abort();
    }
#endif
    return data[i * ncols + j];
}
void Matrix::LimitSize(size_t nRows, size_t nCols)
{
#ifdef gDEBUG
    if(nRows > nrows0 || nCols > ncols0)
    {
        std::cout << "\n Given size exceeds the original size of the matrix" << std::endl;
        abort();
    }
#endif
    nrows = nRows;
    ncols = nCols;
    size = nrows * ncols;
}

/**************************************************************************
   MathLib-Method:
   Task:
   Programing:
   08/2004 WW Implementation
   02/2005 WW Change name
**************************************************************************/
void Matrix::Write(std::ostream& os)
{
    os.setf(std::ios::scientific, std::ios::floatfield);
    os.precision(12);

    for (size_t i = 0; i < nrows; i++)
    {
        os << "| ";
        for (size_t j = 0; j < ncols; j++)
            os << (*this)(i, j) << " ";
        os << "| " << std::endl;
    }
    os << std::endl;
}

/**************************************************************************
   MathLib-Method:
   Task:
   Programing:
   01/2006 WW Implementation
   03/2010 TF write whole matrix in one chunk
**************************************************************************/
void Matrix::Write_BIN(std::fstream& os)
{
    os.write((char*)data, size * sizeof(double));
}
/**************************************************************************
   MathLib-Method:
   Task:
   Programing:
   01/2006 WW Implementation
**************************************************************************/
void Matrix::Read_BIN(std::fstream& is)
{
    for(size_t i = 0; i < size; i++)
        is.read((char*)(&data[i]), sizeof(data[i]));
}

//-----------------------------------------------------
// Symmetrical matrix
SymMatrix::SymMatrix(size_t dim) :
    Matrix(0)
{
    Sym = true;
    nrows = ncols = dim;
    size = (int)nrows * (nrows + 1) / 2;
    data = new double[size];
    nrows0 = ncols0 = dim;
    for(size_t i = 0; i < size; i++)
        data[i] = 0.0;
}

SymMatrix::SymMatrix() : Matrix(0)
{
    Sym = true;
    nrows = 0;
    ncols = 0;
    nrows0 = 0;
    ncols0 = 0;
    size = 0;
    data = 0;
}
SymMatrix::SymMatrix(const SymMatrix& m) : Matrix(0)
{
    Sym = m.Sym;
    nrows = m.nrows;
    ncols = m.ncols;
    nrows0 = m.nrows0;
    ncols0 = m.ncols0;
    size = m.size;
    data = new double[size];
    for(size_t i = 0; i < size; i++)
        data[i] = 0.0;
}

void SymMatrix::resize(size_t dim)
{
    if (size > 0)
    {
        delete[] data;
        data = NULL;
    }

    Sym = true;
    nrows = ncols = dim;
    size = (int) nrows * (nrows + 1) / 2;
    data = new double[size];
    nrows0 = ncols0 = dim;
    for (size_t i = 0; i < size; i++)
        data[i] = 0.0;
}

void SymMatrix::operator = (double a)
{
    for(size_t i = 0; i < size; i++)
        data[i] = a;
}
void SymMatrix::operator *= (double a)
{
    for(size_t i = 0; i < size; i++)
        data[i] *= a;
}
void SymMatrix::operator += (double a)
{
    for(size_t i = 0; i < size; i++)
        data[i] += a;
}

//
void SymMatrix::operator = (const SymMatrix& m)
{
#ifdef gDEBUG
    if(nrows != m.Rows() || ncols != m.Cols())
    {
        std::cout << "\n The sizes of the two matrices are not matched" << endl;
        abort();
    }
#endif
    size_t id = 0;
    for (size_t i = 0; i < nrows; i++)
        for (size_t j = 0; j < ncols; j++)
        {
            if (j > i)
                continue;
            id = i * (i + 1) / 2 + j; // temporary
            data[id] = m(i, j);
        }
}

//
void SymMatrix::operator += (const SymMatrix& m)
{
#ifdef gDEBUG
    if(nrows != m.Rows())
    {
        std::cout << "\n The sizes of the two matrices are not matched" << endl;
        abort();
    }
#endif
    size_t id = 0;
    for (size_t i = 0; i < nrows; i++)
        for (size_t j = 0; j < ncols; j++)
        {
            if (j > i)
                continue;
            id = i * (i + 1) / 2 + j; // temporary
            data[id] += m(i, j);
        }
}

//
void SymMatrix::operator -= (const SymMatrix& m)
{
#ifdef gDEBUG
    if(nrows != m.Rows())                 //Assertion, will be removed
    {
        std::cout << "\n The sizes of the two matrices are not matched" << endl;
        abort();
    }
#endif
    size_t id = 0;
    for (size_t i = 0; i < nrows; i++)
        for (size_t j = 0; j < ncols; j++)
        {
            if (j > i)
                continue;
            id = i * (i + 1) / 2 + j; // temporary
            data[id] -= m(i, j);
        }
}
//
double& SymMatrix::operator() (size_t i, size_t j) const
{
#ifdef gDEBUG
    if(i >= nrows || j >= nrows)
    {
        std::cout << "\n Index exceeds the size of the matrix" << std::endl;
        abort();
    }
#endif

    if(i >= j)
        return data[i * (i + 1) / 2 + j];
    else
        return data[j * (j + 1) / 2 + i];
}

void SymMatrix::LimitSize(size_t dim)
{
#ifdef gDEBUG
    if(dim > nrows0)
    {
        std::cout << "\n Given size exceeds the original size of the matrix" << std::endl;
        abort();
    }
#endif
    nrows = ncols = dim;
    size = nrows * (nrows + 1) / 2;
}

//-----------------------------------------------------
// Diagonal matrix
DiagonalMatrix::DiagonalMatrix(size_t dim) : Matrix(0)
{
    Sym = true;
    nrows = ncols = dim;
    size = dim;
    data = new double[dim];
    nrows0 = ncols0 = dim;
    for(size_t i = 0; i < size; i++)
        data[i] = 0.0;
    dummy_zero = 0.0;
}

DiagonalMatrix::DiagonalMatrix() : Matrix(0)
{
    Sym = true;
    nrows = 0;
    ncols = 0;
    nrows0 = 0;
    ncols0 = 0;
    size = 0;
    data = 0;
    dummy_zero = 0.0;
}
DiagonalMatrix::DiagonalMatrix(const DiagonalMatrix& m) : Matrix(0)
{
    Sym = m.Sym;
    nrows = m.nrows;
    ncols = m.ncols;
    nrows0 = m.nrows0;
    ncols0 = m.ncols0;
    size = m.size;
    data = new double[size];
    for(size_t i = 0; i < size; i++)
        data[i] = 0.0;
    dummy_zero = 0.0;
}

void DiagonalMatrix::resize(size_t dim)
{
    if(size > 0)
    {
        delete [] data;
        data = NULL;
    }

    Sym = true;
    nrows = ncols = dim;
    size = dim;
    data = new double[size];
    nrows0 = ncols0 = dim;
    for(size_t i = 0; i < size; i++)
        data[i] = 0.0;
}

void DiagonalMatrix::operator = (double a)
{
    for(size_t i = 0; i < size; i++)
        data[i] = a;
}
void DiagonalMatrix::operator *= (double a)
{
    for(size_t i = 0; i < size; i++)
        data[i] *= a;
}
void DiagonalMatrix::operator += (double a)
{
    for(size_t i = 0; i < size; i++)
        data[i] += a;
}

//
void DiagonalMatrix::operator = (const DiagonalMatrix& m)
{
#ifdef gDEBUG
    if(nrows != m.Rows() || ncols != m.Cols())
    {
        cout << "\n The sizes of the two matrices are not matched" << endl;
        abort();
    }
#endif
    for(size_t i = 0; i < size; i++)
        data[i] = m(i);
}

//
void DiagonalMatrix::operator += (const DiagonalMatrix& m)
{
#ifdef gDEBUG
    if(nrows != m.Rows())
    {
        cout << "\n The sizes of the two matrices are not matched" << endl;
        abort();
    }
#endif
    for(size_t i = 0; i < size; i++)
        data[i] += m(i);
}

//
void DiagonalMatrix::operator -= (const DiagonalMatrix& m)
{
#ifdef gDEBUG
    if(nrows != m.Rows())                 //Assertion, will be removed
    {
        cout << "\n The sizes of the two matrices are not matched" << endl;
        abort();
    }
#endif
    for(size_t i = 0; i < size; i++)
        data[i] -= m(i);
}
//
double& DiagonalMatrix::operator() (size_t i, size_t j) const
{
#ifdef gDEBUG
    if(i >= nrows || j >= nrows)
    {
        cout << "\n Index exceeds the size of the matrix" << endl;
        abort();
    }
#endif

    if(i == j)
        return data[i];           // temporary
    else
        return dummy_zero;
}

double& DiagonalMatrix::operator() (size_t i) const
{
#ifdef gDEBUG
    if(i >= size)
    {
        cout << "\n Index exceeds the size of the matrix" << endl;
        abort();
    }
#endif

    return data[i];
}

void DiagonalMatrix::LimitSize(size_t dim)
{
#ifdef gDEBUG
    if(dim > nrows0)
    {
        cout << "\n Given size exceeds the original size of the matrix" << endl;
        abort();
    }
#endif
    nrows = ncols = dim;
    size = dim;
}

}
}

// End of class Matrix
//==========================================================================
