/**
 * \file
 * \author Karsten Rink
 * \date   2012-09-27
 * \brief  Implementation of readMeshFromFile function.
 *
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * @file readMeshFromFile.cpp
 * @date 2012-09-27
 * @author Karsten Rink
 */

#include "readMeshFromFile.h"

#ifdef USE_PETSC
#include <petsc.h>
#endif

#include <logog/include/logog.hpp>

#include "BaseLib/FileTools.h"
#include "BaseLib/StringTools.h"

#include "MeshLib/Mesh.h"

#include "FileIO/Legacy/MeshIO.h"
#include "FileIO/VtkIO/VtuInterface.h"

#ifdef USE_MPI
#include "FileIO/MPI_IO/NodePartitionedMeshReader.h"
#include "MeshLib/NodePartitionedMesh.h"
#endif

namespace FileIO
{
MeshLib::Mesh* readMeshFromFile(const std::string &file_name)
{
#ifdef USE_MPI
	NodePartitionedMeshReader read_pmesh(PETSC_COMM_WORLD);
	return read_pmesh.read(file_name);
#else
	if (BaseLib::hasFileExtension("msh", file_name))
	{
		Legacy::MeshIO meshIO;
		return meshIO.loadMeshFromFile(file_name);
	}

	if (BaseLib::hasFileExtension("vtu", file_name))
		return VtuInterface::readVTUFile(file_name);

	ERR("readMeshFromFile(): Unknown mesh file format in file %s.", file_name.c_str());
	return nullptr;
#endif
}
} // end namespace FileIO
