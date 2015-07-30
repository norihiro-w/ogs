/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file PVDOutput.cpp
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#include "PVDOutput.h"

#include <string>
#include <vector>

#include <logog/include/logog.hpp>

#include "BaseLib/MPITools.h"

#include "FileIO/OutputIO/PVDWriter.h"
#include "FileIO/OutputIO/VtuWriter.h"
#include "FileIO/OutputIO/PVtuWriter.h"


void PVDOutput::write(  const NumLib::TimeStep &current_time, 
                        BaseLib::OrderedMap<std::string, OutputVariableInfo> &data)
{
    BaseLib::MPIEnvironment mpi;

    // prepare vtu data
    std::vector<VtuWriter::PointData> node_values;
    std::vector<VtuWriter::CellData> ele_values;
    
    for (BaseLib::OrderedMap<std::string, OutputVariableInfo>::iterator itr = data.begin(); itr!=data.end(); ++itr) {
        // pick up variables for this mesh and specified by users
        if (itr->second.mesh_id != getMesh()->getID()) continue;
        if ((itr->second.object_type == OutputVariableInfo::Node && hasNodalVariable(itr->first))
            || (itr->second.object_type == OutputVariableInfo::Element && hasElementalVariable(itr->first))) {
            OutputVariableInfo &var = itr->second;
            VtuWriter::AttributeInfo attr(var.name, var.nr_of_components, var.value);
            switch (var.data_type) {
                case OutputVariableInfo::Char: attr.data_type = VtuWriter::Char; break;
                case OutputVariableInfo::Int: attr.data_type = VtuWriter::Int; break;
                case OutputVariableInfo::Real: attr.data_type = VtuWriter::Real; break;
            }
            if (var.object_type==OutputVariableInfo::Node) {
                node_values.push_back(VtuWriter::PointData(var.name, attr));
            } else if (var.object_type==OutputVariableInfo::Element) {
                ele_values.push_back(VtuWriter::CellData(var.name, attr));
            }
        }
    }
    
    if (node_values.size() == 0 && ele_values.size() == 0) {
        WARN("***Warning: Asked to write results but specified data not found. Skip this output. Please check consistency of variable names in process and output settings.");
        return;
    }

    // set base file name for output
    const std::string vtk_file_basename = getOutputBaseName() + "_" + std::to_string(current_time.steps());
    std::string vtk_file_path = getOutputDir();
    if (vtk_file_path.length()>0) vtk_file_path += "/";

    // write VTU file
    std::string vtu_file_name_base = vtk_file_basename;
    if (mpi.isParallel())
        vtu_file_name_base += "_part" + std::to_string(mpi.rank());
    vtu_file_name_base += ".vtu";
    std::string vtu_file_name = vtk_file_path + vtu_file_name_base;
    INFOa("Writing results...: %s", vtu_file_name.c_str());
    VtuWriter vtuWriter(false);
    vtuWriter.write(vtu_file_name, *getMesh(), node_values, ele_values, true);
    
    if (mpi.rank() == 0)
    {
        std::string vtkfile = vtu_file_name;
        if (mpi.isParallel())
        {
            vtkfile = vtk_file_basename + ".pvtu";
            std::string pvtu_file_name = vtk_file_path + vtk_file_basename + ".pvtu";
            PVtuWriter pvtu;
            pvtu.write(pvtu_file_name, *getMesh(), node_values, ele_values, true);
        }

        // update PVD file
        if (_pvd==NULL) {
            _pvd = new PVDWriter();
            std::string pvd_name = getOutputDir();
            if (pvd_name.length()>0) pvd_name += "/";
            pvd_name += getOutputBaseName() + ".pvd";
            //std::string pvd_name = getOutputDir()+"\\"+getOutputBaseName() + ".pvd";
            _pvd->initialize(pvd_name);
        }

        _pvd->update(current_time.current(), vtkfile);
    }
}
