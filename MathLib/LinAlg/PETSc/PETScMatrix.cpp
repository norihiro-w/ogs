/*!
   \file  PETScMatrix.cpp
   \brief Definition of member functions of class PETScMatrix, which provides an interface to
          PETSc matrix routines.

   \author Wenqing Wang
   \date Nov 2011 - Sep 2013

   \copyright
   Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
               Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license
*/

#include "PETScMatrix.h"

#include <logog/include/logog.hpp>

#include <petscmat.h>

namespace MathLib
{

PETScMatrix::PETScMatrix (const PetscInt nrows)
    : _external_data(false)
{
    PetscInt m = PETSC_DECIDE, n = PETSC_DECIDE, M = nrows, N = PETSC_DECIDE;

    MatCreate(PETSC_COMM_WORLD, &_A);
    MatSetSizes(_A, m, n, M, N);

//#ifdef USE_MPI
//    MatSetType(_A, MATMPIAIJ);
//#else
//    MatSetType(_A, MATSEQAIJ);
//#endif
    MatSetFromOptions(_A);

    MatGetOwnershipRange(_A, &_start_rank, &_end_rank);
    MatGetSize(_A, &_nrows,  &_ncols);
    MatGetLocalSize(_A, &_n_loc_rows, &_n_loc_cols);
}

PETScMatrix::PETScMatrix (const PetscInt nrows, const MatrixOption& mat_opt)
    : _external_data(false)
{
    PetscInt m = PETSC_DECIDE, n = PETSC_DECIDE, M = nrows, N = PETSC_DECIDE;
    if (!mat_opt.is_global_size)
    {
//        int local_rows = m;
//        int all_rows = 0;
//        MPI_Allreduce(&local_rows, &all_rows, 1, MPI_INT, MPI_SUM, PETSC_COMM_WORLD);
//        //N = M; // assume square
        m = nrows;
        M = PETSC_DECIDE;
        n = m; // n should correspond to the local size of solution and rhs vectors
        N = PETSC_DECIDE; //assume square matrix
    }

    MatCreate(PETSC_COMM_WORLD, &_A);
    MatSetSizes(_A, m, n, M, N);

//#ifdef USE_MPI
//    MatSetType(_A, MATMPIAIJ);
//#else
//    MatSetType(_A, MATSEQAIJ);
//#endif
    MatSetFromOptions(_A);

    PetscInt d_nz = (mat_opt.d_nz == 0) ? PETSC_DECIDE : mat_opt.d_nz;
    PetscInt o_nz = (mat_opt.o_nz == 0) ? PETSC_DECIDE : mat_opt.o_nz;
    MatSeqAIJSetPreallocation(_A, d_nz, PETSC_NULL);
    MatMPIAIJSetPreallocation(_A, d_nz, PETSC_NULL, o_nz, PETSC_NULL);
    // If pre-allocation does not work one can use MatSetUp(_A), which is much
    // slower.

    MatGetOwnershipRange(_A, &_start_rank, &_end_rank);
    MatGetSize(_A, &_nrows,  &_ncols);
    MatGetLocalSize(_A, &_n_loc_rows, &_n_loc_cols);

    INFO("creating a matrix: rows=%d, cols=%d, lc rows=%d, lc cols=%d, start=%d, end=%d", _nrows, _ncols, _n_loc_rows, _n_loc_cols, _start_rank, _end_rank);
}

//PETScMatrix::PETScMatrix (const PetscInt nrows, const PetscInt ncols, const MatrixOption* mat_opt)
//    :_nrows(nrows), _ncols(ncols),  _n_loc_rows(PETSC_DECIDE),
//     _n_loc_cols(mat_opt.n_local_cols), _external_data(false)
//{
//    if(!mat_opt.is_global_size)
//    {
//        _nrows = PETSC_DECIDE;
//        _ncols = PETSC_DECIDE;
//        _n_loc_rows = nrows;
//        _n_loc_cols = ncols;
//    }
//
//    create(mat_opt.d_nz, mat_opt.o_nz);
//}

PETScMatrix::PETScMatrix (Mat &mat)
: _A(mat), _external_data(true)
{
    MatGetSize(_A, &_nrows, &_ncols);
    MatGetLocalSize(_A, &_n_loc_rows, &_n_loc_cols);
    MatGetOwnershipRange(_A, &_start_rank, &_end_rank);
}

void PETScMatrix::setRowsColumnsZero(std::vector<std::size_t> const& row_pos)
{
    // Each rank (compute core) processes only the rows that belong to the rank itself.
    const PetscScalar one = 1.0;
    const PetscInt nrows = static_cast<PetscInt> (row_pos.size());

    if(nrows>0)
        MatZeroRows(_A, nrows, (PetscInt*)&row_pos[0], one, PETSC_NULL, PETSC_NULL);
    else
        MatZeroRows(_A, 0, PETSC_NULL, one, PETSC_NULL, PETSC_NULL);
}

void PETScMatrix::viewer(const std::string &file_name, const PetscViewerFormat vw_format) const
{
    PetscViewer viewer;
    PetscViewerASCIIOpen(PETSC_COMM_WORLD, file_name.c_str(), &viewer);
    PetscViewerPushFormat(viewer, vw_format);

    finalizeAssembly();

    PetscObjectSetName((PetscObject)_A,"Stiffness_matrix");
    MatView(_A,viewer);

    PetscViewerDestroy(&viewer);

// This preprocessor is only for debugging, e.g. dump the matrix and exit the program.
//#define EXIT_TEST
#ifdef EXIT_TEST
    MatDestroy(&_A);
    PetscFinalize();
    exit(0);
#endif

}

//void PETScMatrix::create(const PetscInt d_nz, const PetscInt o_nz)
//{
//    MatCreate(PETSC_COMM_WORLD, &_A);
//    MatSetSizes(_A, _n_loc_rows, _n_loc_cols, _nrows, _ncols);
//
//#ifdef USE_MPI
//    MatSetType(_A, MATMPIAIJ);
//#else
//    MatSetType(_A, MATSEQAIJ);
//#endif
//    MatSetFromOptions(_A);
//
//    MatSeqAIJSetPreallocation(_A, d_nz, PETSC_NULL);
//    MatMPIAIJSetPreallocation(_A, d_nz, PETSC_NULL, o_nz, PETSC_NULL);
//    // If pre-allocation does not work one can use MatSetUp(_A), which is much
//    // slower.
//
//    MatGetOwnershipRange(_A, &_start_rank, &_end_rank);
//    MatGetSize(_A, &_nrows,  &_ncols);
//    MatGetLocalSize(_A, &_n_loc_rows, &_n_loc_cols);
//}

bool finalizeMatrixAssembly(PETScMatrix &mat, const MatAssemblyType asm_type)
{
    mat.finalizeAssembly(asm_type);
    return true;
}

} //end of namespace

