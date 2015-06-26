/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file PVDWriter.cpp
 *
 * Created on 2012-07-17 by Norihiro Watanabe
 */


#include "PVDWriter.h"

#include <fstream>
#include "logog/include/logog.hpp"

const std::string INDEX_STR = "  ";

bool PVDWriter::initialize(const std::string &pvd_file)
{
    //PVD
    this->_pvd_data.vec_dataset.clear();
    this->_pvd_file_name = pvd_file;

    return true;
}

bool PVDWriter::update(double time, const std::string &vtk_file)
{
    PVDData::VTK_Info vtk_info;
    vtk_info.timestep = time;
    vtk_info.vtk_file = vtk_file;
    this->_pvd_data.vec_dataset.push_back(vtk_info);
    
    std::fstream fin(this->_pvd_file_name.c_str(), std::ios::out);
    if (!fin.good()) {
        WARN("***Warning: Cannot open the PVD file, %s", this->_pvd_file_name.c_str());
        return false;
    }

    PVDWriter::writeHeader(fin);
    for (int i = 0; i < (int)_pvd_data.vec_dataset.size(); i++)
        PVDWriter::writeDataset(fin, _pvd_data.vec_dataset[i].timestep, _pvd_data.vec_dataset[i].vtk_file);
    PVDWriter::writeEnd(fin);

    fin.close();

    return true;
    
}

bool PVDWriter::writeHeader(std::fstream &fin)
{
    fin << "<?xml version=\"1.0\"?>" << std::endl;
    fin <<
    "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">"
        << std::endl;
    fin << INDEX_STR << "<Collection>" << std::endl;
    return true;
}

bool PVDWriter::writeEnd(std::fstream &fin)
{
    fin << INDEX_STR << "</Collection>" << std::endl;
    fin << "</VTKFile>" << std::endl;
    return true;
}

bool PVDWriter::writeDataset(std::fstream &fin, double timestep, const std::string &vtkfile)
{
    fin.setf(std::ios::scientific,std::ios::floatfield);
    fin.precision(12);
    fin << INDEX_STR << INDEX_STR << "<DataSet timestep=\"" << timestep <<
    "\" group=\"\" part=\"0\" file=\"" << vtkfile << "\"/>" << std::endl;
    return true;
}

