/**
 * \file  PETScVector.cpp
 * \brief Definition of member functions of class PETScVector,
 *        which provides an interface to PETSc vector routines.
 *
 *   Note: the return message of PETSc routines is ommited in
 *         the source code. If it is really needed, it can be activated by
 *         adding a PetscErrorCode type variable before each PETSc fucntion
 *
 * \author Wenqing Wang
 * \date Nov 2011 - Sep 2013
 *
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */


#include "PETScVector.h"

#include <logog/include/logog.hpp>

namespace MathLib
{

PETScVector::PETScVector(const PetscInt vec_size, const bool is_sequential, const std::vector<std::size_t>* vec_ghosts)
: _lv(nullptr), _external_data(false)
{
    if (is_sequential) {
        VecCreateSeq(PETSC_COMM_SELF, vec_size, &_v);
    } else if (vec_ghosts == nullptr || vec_ghosts->empty()) {
        VecCreateMPI(PETSC_COMM_WORLD, vec_size, PETSC_DECIDE, &_v);
    } else {
        PetscInt nghost = vec_ghosts->size();
        PetscInt* ifrom = (PetscInt*)&(*vec_ghosts)[0];
        VecCreateGhost(PETSC_COMM_WORLD, vec_size, PETSC_DECIDE, nghost, ifrom, &_v);
        VecGhostGetLocalForm(_v, &_lv);
        for (auto i : *vec_ghosts)
            _vec_ghosts_gid.push_back(i);
    }

    VecGetOwnershipRange(_v, &_start_rank, &_end_rank);
    VecGetLocalSize(_v, &_size_loc);
    VecGetSize(_v, &_size);

    INFO("creating a vector: size=%d, local size=%d, start=%d, end=%d", _size, _size_loc, _start_rank, _end_rank);

    assert(sizeof(PetscInt) == sizeof(std::size_t)); // to cast size_t* to PetscInt*
}

PETScVector::PETScVector(const PETScVector &existing_vec, const bool deep_copy)
: _external_data(false), _vec_ghosts_gid(existing_vec. _vec_ghosts_gid)
{
    VecDuplicate(existing_vec._v, &_v);
    VecGhostGetLocalForm(_v, &_lv);

    VecGetOwnershipRange(_v, &_start_rank,&_end_rank);
    VecGetLocalSize(_v, &_size_loc);
    VecGetSize(_v, &_size);

    // Copy values
    if(deep_copy)
    {
        finalizeAssembly();
        VecCopy(existing_vec._v, _v);
    }
}


PETScVector::PETScVector(Vec &vec)
: _v(vec), _lv(nullptr), _external_data(true)
{
    VecGetOwnershipRange(_v, &_start_rank,&_end_rank);
    VecGetLocalSize(_v, &_size_loc);
    VecGetSize(_v, &_size);
}

PETScVector::~PETScVector()
{
    if (_external_data) return;
    if (_lv)
        VecGhostRestoreLocalForm(_v,& _lv);
    VecDestroy(&_v);
}

void PETScVector::finalizeAssembly()
{
    VecAssemblyBegin(_v);
    VecAssemblyEnd(_v);

    if (_lv!=nullptr) {
        VecGhostUpdateBegin(_v, INSERT_VALUES,SCATTER_FORWARD);
        VecGhostUpdateEnd(_v, INSERT_VALUES,SCATTER_FORWARD);
    }
}

IVector& PETScVector::operator = (const IVector &v_in_)
{
    const PETScVector &v_in(static_cast<const PETScVector &>(v_in_));
    const_cast<PETScVector&>(v_in).finalizeAssembly();
    VecCopy(v_in._v, _v);
    finalizeAssembly();
    return *this;
}

void PETScVector::gatherLocalVectors( PetscScalar local_array[],
                                      PetscScalar global_array[])
{
    // Collect vectors from processors.
    int size_rank;
    MPI_Comm_size(PETSC_COMM_WORLD, &size_rank);

    // number of elements to be sent for each rank
    std::vector<int>  i_cnt(size_rank);

    MPI_Allgather(&_size_loc, 1, MPI_INT, &i_cnt[0], 1, MPI_INT, PETSC_COMM_WORLD);

    // collect local array
    int offset = 0;
    // offset in the receive vector of the data from each rank
    std::vector<int>  i_disp(size_rank);
    for(int i=0; i<size_rank; i++)
    {
        i_disp[i] = offset;
        offset += i_cnt[i];
    }

    MPI_Allgatherv(local_array, _size_loc, MPI_DOUBLE,
                   global_array, &i_cnt[0], &i_disp[0], MPI_DOUBLE, PETSC_COMM_WORLD);

}

void PETScVector::getGlobalVector(PetscScalar u[])
{

#ifdef TEST_MEM_PETSC
    PetscLogDouble mem1, mem2;
    PetscMemoryGetCurrentUsage(&mem1);
#endif

    PetscScalar *xp = nullptr;
    VecGetArray(_v, &xp);

    gatherLocalVectors(xp, u);

    //This following line may be needed late on
    //  for a communication load balance:
    //MPI_Barrier(PETSC_COMM_WORLD);

    VecRestoreArray(_v, &xp);

    //TEST
#ifdef TEST_MEM_PETSC
    PetscMemoryGetCurrentUsage(&mem2);
    PetscPrintf(PETSC_COMM_WORLD, "### Memory usage by Updating. Before :%f After:%f Increase:%d\n", mem1, mem2, (int)(mem2 - mem1));
#endif
}

PetscScalar PETScVector::getNorm(MathLib::VecNormType nmtype) const
{
    NormType petsc_norm = NORM_1;
    switch(nmtype)
    {
        case MathLib::VecNormType::NORM1:
            petsc_norm = NORM_1;
            break;
        case MathLib::VecNormType::NORM2:
            petsc_norm = NORM_2;
            break;
        case MathLib::VecNormType::INFINITY_N:
            petsc_norm = NORM_INFINITY;
            break;
        default:
            break;
    }

    PetscScalar norm = 0.;
    VecNorm(_v, petsc_norm, &norm);
    return norm;
}

void PETScVector::viewer(const std::string &file_name, const PetscViewerFormat vw_format) const
{
    PetscViewer viewer;
    PetscViewerASCIIOpen(PETSC_COMM_WORLD, file_name.c_str(), &viewer);
    PetscViewerPushFormat(viewer, vw_format);

    PetscObjectSetName((PetscObject)_v, file_name.c_str());
    VecView(_v, viewer);

#define  nEXIT_TEST
#ifdef EXIT_TEST
    VecDestroy(&_v);
    PetscFinalize();
    exit(0);
#endif

}

void finalizeVectorAssembly(PETScVector &vec)
{
    vec.finalizeAssembly();
}

} //end of namespace
