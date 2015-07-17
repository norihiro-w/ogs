/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 */

#include "PVtuWriter.h"

#include <fstream>

#include "logog/include/logog.hpp"

#include "BaseLib/SystemTools.h"
#include "BaseLib/FileTools.h"
#include "BaseLib/MPITools.h"

#include "MeshLib/CoordinateSystem.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"
#include "MeshLib/Elements/Element.h"

using namespace MeshLib;

PVtuWriter::PVtuWriter()
: _isInitialized(false), _isLittleEndian(true),
  type_UChar(VtuWriter::Int8), type_Int(VtuWriter::Int32),
  type_UInt(VtuWriter::Int8), type_Long(VtuWriter::Int8),
  type_Double(VtuWriter::Int8),  SIZE_OF_BLOCK_LENGTH_TAG(0)
{
    this->initialize();
}

bool PVtuWriter::write(const std::string &vtkfile,
                                    MeshLib::Mesh &msh, std::vector<PointData> &nodal_values,
                                    std::vector<CellData> &elemental_values, bool outGroupID)
{
    //-------------------------------------------------------------------------
    //# Setup file stream
    //-------------------------------------------------------------------------
    std::fstream fin;
    fin.open(vtkfile.c_str(), std::ios::out);

    if (!fin.good())
    {
        WARN("***Warning: Cannot open the output file, %s", vtkfile.c_str());
        return false;
    }

    fin.setf(std::ios::scientific,std::ios::floatfield);
    fin.precision(12);

    for (size_t i=0; i<nodal_values.size(); i++) {
        VtuWriter::AttributeInfo &attr = nodal_values[i].second;
        switch (attr.data_type) {
        case VtuWriter::Char: attr.vtk_data_type = type_UChar; break;
        case VtuWriter::Int: attr.vtk_data_type = type_Long; break;
        case VtuWriter::Real: attr.vtk_data_type = type_Double; break;
        }
    }
    for (size_t i=0; i<elemental_values.size(); i++) {
        VtuWriter::AttributeInfo &attr = elemental_values[i].second;
        switch (attr.data_type) {
        case VtuWriter::Char: attr.vtk_data_type = type_UChar; break;
        case VtuWriter::Int: attr.vtk_data_type = type_Long; break;
        case VtuWriter::Real: attr.vtk_data_type = type_Double; break;
        }
    }

    //-------------------------------------------------------------------------
    //# Output
    //-------------------------------------------------------------------------
    long offset = 0;

    const std::string str_format = "ascii";

    //# Header
    fin << "<?xml version=\"1.0\"?>" << std::endl;
    fin << "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\"";
    if (_isLittleEndian)
        fin << " byte_order=\"LittleEndian\"";
    else
        fin << " byte_order=\"BigEndian\"";
    fin << ">" << std::endl;
    //  fin << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">"  << endl;

    //# Unstructured Grid information
    fin << "  <PUnstructuredGrid GhostLevel=\"1\">\n";
    fin << "    <PPoints>\n";
    writeDataArrayHeader(fin, type_Double, "", 3, str_format, offset);
    fin << "    </PPoints>\n";

    if (!nodal_values.empty())
    {
        fin << "    <PPointData Scalars=\"" << nodal_values[0].first << "\" >\n";
        writePointData(fin, true, nodal_values, msh, offset);
        fin << "    </PPointData>\n";
    }

    if (elemental_values.empty())
        fin << "    <PCellData Scalars=\"" << "MatGroup" << "\" >\n";
    else
        fin << "    <PCellData Scalars=\"" << elemental_values[0].first << "\" >\n";
//    writeDataArrayHeader(fin, type_Int, "Domain", 1, str_format, offset);
    writeCellData(fin, true, elemental_values, msh, offset);
    writeDataArrayHeader(fin, type_Int, "MatGroup", 1, str_format, offset);
    fin << "    </PCellData>\n";

    // sub files
    BaseLib::MPIEnvironment mpi;
    const std::string path_base = BaseLib::extractPath(vtkfile) + BaseLib::extractBaseNameWithoutExtension(vtkfile);
    for (int i=0; i<mpi.size(); i++) {
        std::string sub_vtu_file_name = path_base + "_part" + std::to_string(i) + ".vtu";
        fin << "    <Piece Source=\"" << sub_vtu_file_name << "\" />\n";
    }
    // closing
    fin << "  </PUnstructuredGrid>" << "\n";
    fin << "</VTKFile>" << "\n";
    fin.close();

    return true;
}

bool PVtuWriter::writePointData(std::fstream &fin,
                           bool /*output_data*/,
                           std::vector<PointData> &nod_values,
                           Mesh &msh,
                           long &offset)
{
    //Nodal values
    for (size_t i=0; i<nod_values.size(); i++)
    {
        const std::string &var_name = nod_values[i].first;
        auto &pt_data = nod_values[i].second;
        writeDataArrayHeader(fin, pt_data.vtk_data_type, var_name, pt_data.nr_of_components, getFormatType(), offset);
    }

    return true;
}

bool PVtuWriter::writeCellData(std::fstream &fin,
                             bool /*output_data*/,
                             std::vector<CellData> &ele_values,
                             Mesh &msh,
                             long &offset)
{
    //Element values
    for (int i = 0; i < (int) ele_values.size(); i++)
    {
        auto &pt_data = ele_values[i].second;
        writeDataArrayHeader(fin, pt_data.vtk_data_type, pt_data.attribute_name, pt_data.nr_of_components, getFormatType(), offset);
    }

    return true;
}


void PVtuWriter::initialize()
{
    //if (useBinary) {
    //======================================================================
    //# Set machine dependent stuff
    //Data type
    if (sizeof(unsigned char) == 1)
        type_UChar = VtuWriter::UInt8;
    else if (sizeof(unsigned char) == 2)
        type_UChar = VtuWriter::UInt16;
    if (sizeof(int) == 4)
        type_Int = VtuWriter::Int32;
    else if (sizeof(int) == 8)
        type_Int = VtuWriter::Int64;
    if (sizeof(unsigned int) == 4)
        type_UInt = VtuWriter::UInt32;
    else if (sizeof(unsigned int) == 8)
        type_UInt = VtuWriter::UInt64;
    if (sizeof(long) == 4)
        type_Long = VtuWriter::Int32;
    else if (sizeof(long) == 8)
        type_Long = VtuWriter::Int64;
    if (sizeof(double) == 4)
        type_Double = VtuWriter::Float32;
    else if (sizeof(double) == 8)
        type_Double = VtuWriter::Float64;
    //
    SIZE_OF_BLOCK_LENGTH_TAG = sizeof(unsigned int);
    //Endian(byte order)
    _isLittleEndian = BaseLib::IsLittleEndian();
    //}

    this->_isInitialized = true;
}

size_t PVtuWriter::getSizeOfVtkDataType(VtuWriter::DataType data_type) const
{
    size_t n = 0;
    switch (data_type) {
    case VtuWriter::Char: n = sizeof(unsigned char); break;
    case VtuWriter::Int: n = sizeof(long); break;
    case VtuWriter::Real: n = sizeof(double); break;
    }

    return n;
}

bool PVtuWriter::writeDataArrayHeader(std::fstream &fin,
                                VtuWriter::VTK_XML_DATA_TYPE data_type,
                                const std::string &str_name,
                                int nr_components,
                                const std::string &str_format,
                                long offset)
{
    std::string str_data_type;
    switch (data_type)
    {
    case VtuWriter::Int8: str_data_type = "Int8";
        break;
    case VtuWriter::UInt8: str_data_type = "UInt8";
        break;
    case VtuWriter::Int16: str_data_type = "Int16";
        break;
    case VtuWriter::UInt16: str_data_type = "UInt16";
        break;
    case VtuWriter::Int32: str_data_type = "Int32";
        break;
    case VtuWriter::UInt32: str_data_type = "UInt32";
        break;
    case VtuWriter::Int64: str_data_type = "Int64";
        break;
    case VtuWriter::UInt64: str_data_type = "UInt64";
        break;
    case VtuWriter::Float32: str_data_type = "Float32";
        break;
    case VtuWriter::Float64: str_data_type = "Float64";
        break;
    }
    fin << "        <PDataArray type=\"" << str_data_type << "\"";
    if (str_name != "")
        fin << " Name=\"" << str_name << "\"";
    if (nr_components > 1)
        fin << " NumberOfComponents=\"" << nr_components << "\"";
    fin << " format=\"" << str_format << "\"";
    fin << " />" << std::endl;

    return true;
}


