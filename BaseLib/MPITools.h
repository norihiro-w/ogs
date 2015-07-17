/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef MPITOOLS_H
#define MPITOOLS_H

#ifdef USE_MPI
#include <mpi.h>
#endif

namespace BaseLib
{

class MPIEnvironment
{
public:
#ifdef USE_MPI
    MPIEnvironment(MPI_Comm comm = MPI_COMM_WORLD)
    : _comm(comm)
    {
        MPI_Comm_rank(comm, &_rank);
        MPI_Comm_size(comm, &_size);
    }
#else
    typedef int MPI_Comm;
    MPIEnvironment(MPI_Comm comm = 0)
    : _comm(comm), _rank(0), _size(1) {}
#endif

    MPI_Comm communicator() const { return _comm; }
    int rank() const { return _rank; }
    int size() const { return _size; }
    bool root() const { return rank()==0; }
    bool isParallel() const { return size() > 1; }
    void barrier() const
    {
#ifdef USE_MPI
        MPI_Barrier(_comm);
#endif
    }
private:
    MPI_Comm _comm;
    int _rank;
    int _size;
};

} // BaseLib

#endif
