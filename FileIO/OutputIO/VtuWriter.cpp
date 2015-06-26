/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file VtuWriter.cpp
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#include "VtuWriter.h"

#include <fstream>

#include "logog/include/logog.hpp"

#include "BaseLib/SystemTools.h"
#include "BaseLib/FileTools.h"

#include "MeshLib/CoordinateSystem.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"
#include "MeshLib/Elements/Element.h"

using namespace MeshLib;

VtuWriter::VtuWriter(bool binary_mode)
: _useBinary(binary_mode), _isInitialized(false), _isLittleEndian(true),
  type_UChar(Int8), type_Int(Int32),
  type_UInt(Int8), type_Long(Int8),
  type_Double(Int8),  SIZE_OF_BLOCK_LENGTH_TAG(0)
{
    this->initialize();
}

bool VtuWriter::write(const std::string &vtkfile,
                                    MeshLib::Mesh &msh, std::vector<PointData> &nodal_values, std::vector<CellData> &elemental_values, bool outGroupID)
{
    //-------------------------------------------------------------------------
    //# Setup file stream
    //-------------------------------------------------------------------------
    std::fstream fin;
    if (_useBinary)
        fin.open(vtkfile.c_str(), std::ios::out | std::ios::binary);
    else
        fin.open(vtkfile.c_str(), std::ios::out);

    if (!fin.good())
    {
        WARN("***Warning: Cannot open the output file, %s", vtkfile.c_str());
        return false;
    }

    if (!_useBinary)
    {
        fin.setf(std::ios::scientific,std::ios::floatfield);
        fin.precision(12);
    }

    for (size_t i=0; i<nodal_values.size(); i++) {
        AttributeInfo &attr = nodal_values[i].second;
        switch (attr.data_type) {
        case VtuWriter::Char: attr.vtk_data_type = type_UChar; break;
        case VtuWriter::Int: attr.vtk_data_type = type_Long; break;
        case VtuWriter::Real: attr.vtk_data_type = type_Double; break;
        }
    }
    for (size_t i=0; i<elemental_values.size(); i++) {
        AttributeInfo &attr = elemental_values[i].second;
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

    std::string str_format;
    if (!_useBinary)
        str_format = "ascii";
    else
        str_format = "appended";
    bool data_out = !_useBinary;

    //# Header
    fin << "<?xml version=\"1.0\"?>" << std::endl;
    fin << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\"";
    if (!_useBinary || _isLittleEndian)
        fin << " byte_order=\"LittleEndian\"";
    else
        fin << " byte_order=\"BigEndian\"";
    fin << ">" << std::endl;
    //  fin << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">"  << endl;

    //# Unstructured Grid information
    fin << "  <UnstructuredGrid>" << std::endl;
    fin << "    <Piece NumberOfPoints=\"" << msh.getNNodes() <<
    "\" NumberOfCells=\"" << msh.getNElements() << "\">" << std::endl;
    //....................................................................
    // Nodes
    fin << "      <Points>" << std::endl;
    writeDataArrayHeader(fin, type_Double, "", 3, str_format, offset);
    writePoints(fin, data_out, msh, offset);
    writeDataArrayFooter(fin);
    fin << "      </Points>" << std::endl;
    //....................................................................
    // Elements
    fin << "      <Cells>" << std::endl;
    //connectivity
    writeDataArrayHeader(fin, type_Long, "connectivity", 0, str_format, offset);
    long sum_ele_components = 0;
    writeCellConnectivity(fin, data_out, msh, offset, sum_ele_components);
    writeDataArrayFooter(fin);
    //offset
    writeDataArrayHeader(fin, type_Long, "offsets", 0, str_format, offset);
    writeCellOffset(fin, data_out, msh, offset);
    writeDataArrayFooter(fin);
    //type
    writeDataArrayHeader(fin, type_UChar, "types", 0, str_format, offset);
    writeCellType(fin, data_out, msh, offset);
    writeDataArrayFooter(fin);
    fin << "      </Cells>" << std::endl;
    //....................................................................
    // Nodal values
    if (nodal_values.size() > 0)
        fin << "      <PointData Scalars=\"" << nodal_values[0].first << "\">" << std::endl;
    else
        fin << "      <PointData Scalars=\"scalars\">" << std::endl;
    writePointData(fin, data_out, nodal_values, msh, offset);
    fin << "      </PointData>" << std::endl;
    //....................................................................
    // Element values
    fin << "      <CellData>" << std::endl;
    writeCellData(fin, data_out, elemental_values, msh, offset);
    if (outGroupID)
        writeElementGroupID(fin, data_out, msh, offset);
    fin << "      </CellData>" << std::endl;
    fin << "    </Piece>" << std::endl;
    fin << "  </UnstructuredGrid>" << std::endl;

    //======================================================================
    // Raw data (for binary mode)
    if (_useBinary)
    {
        fin << "  <AppendedData encoding=\"raw\">" << std::endl;
        fin << "    _";

        //Node
        this->writePoints(fin, true, msh, offset);
        //Element
        //connectivity
        this->writeCellConnectivity(fin, true, msh, offset, sum_ele_components);
        //offset
        this->writeCellOffset(fin, true, msh, offset);
        //type
        this->writeCellType(fin, true, msh, offset);
        // Nodal values
        this->writePointData(fin, true, nodal_values, msh, offset);
        // Elemental values
        this->writeCellData(fin, true, elemental_values, msh, offset);

        fin << std::endl;
        fin << "  </AppendedData>" << std::endl;
    }

    fin << "</VTKFile>" << std::endl;
    fin.close();

    return true;
}

bool VtuWriter::writePointData(std::fstream &fin,
                           bool output_data,
                           std::vector<PointData> &nod_values,
						   Mesh &msh,
                           long &offset)
{
    //Nodal values
    for (size_t i=0; i<nod_values.size(); i++)
    {
        const std::string &var_name = nod_values[i].first;
        AttributeInfo &pt_data = nod_values[i].second;
        unsigned int bin_data_length = getSizeOfVtkDataType(pt_data.data_type)* pt_data.nr_of_components* msh.getNNodes();

        if (!_useBinary || !output_data)
            writeDataArrayHeader(fin, pt_data.vtk_data_type, var_name, pt_data.nr_of_components, getFormatType(), offset);

        if (output_data)
        {
            if (!_useBinary) {
                fin << "          ";
            } else {
                BaseLib::writeValueBinary<unsigned int> (fin, bin_data_length);
            }
            MathLib::LocalMatrix v;
            for (size_t j = 0; j < msh.getNNodes(); j++) {
                pt_data.f->eval(NumLib::TXPosition(NumLib::TXPosition::Node, j, msh.getNode(j)->getCoords()), v);
                if (!_useBinary) {
                    if (v.array().size()>=(int)pt_data.nr_of_components) {
                        for (size_t k=0; k<pt_data.nr_of_components; k++)
                            fin << v.array()(k) << " ";
                    } else {
                        for (size_t k=0; k<pt_data.nr_of_components; k++)
                            fin << 0 << " ";
                    }
                } else {
                    for (size_t k=0; k<pt_data.nr_of_components; k++)
                        BaseLib::writeValueBinary(fin, v.array()(k));
                }
            }
            if (!_useBinary) {
                fin << std::endl;
            }
        }
        else
        {
            offset += bin_data_length + SIZE_OF_BLOCK_LENGTH_TAG;
        }

        if (!_useBinary || !output_data)
            writeDataArrayFooter(fin);
    }


    return true;
}

bool VtuWriter::writeCellData(std::fstream &fin,
                             bool output_data,
                             std::vector<CellData> &ele_values,
							 Mesh &msh,
                             long &offset)
{
    //Element values
    for (int i = 0; i < (int) ele_values.size(); i++)
    {
        AttributeInfo &pt_data = ele_values[i].second;
        unsigned int bin_data_length = getSizeOfVtkDataType(pt_data.data_type)* pt_data.nr_of_components* msh.getNElements();

        if (!_useBinary || !output_data)
            writeDataArrayHeader(fin, pt_data.vtk_data_type, pt_data.attribute_name, pt_data.nr_of_components, getFormatType(), offset);

        if (output_data)
        {
            if (!_useBinary)
            {
                fin << "          ";
                MathLib::LocalMatrix v;
                for (size_t j=0; j<msh.getNElements(); j++) {
                    pt_data.f->eval(NumLib::TXPosition(NumLib::TXPosition::Element, j), v);
                    const size_t n_dummy = pt_data.nr_of_components - static_cast<size_t>(v.array().size());
                    for (int k=0; k<v.array().size(); k++)
                        fin << v.array()(k) << " ";
                    for (size_t k=0; k<n_dummy; k++)
                        fin << 0 << " ";
                }
                fin << std::endl;
            }
            else
            {
                BaseLib::writeValueBinary<unsigned int> (fin, bin_data_length);
                MathLib::LocalMatrix v;
                for (long j = 0; j < (long) msh.getNElements(); j++) {
                    pt_data.f->eval(NumLib::TXPosition(NumLib::TXPosition::Element, j), v);
                    const size_t n_dummy = pt_data.nr_of_components - static_cast<size_t>(v.array().size());
                    for (int k=0; k<v.array().size(); k++)
                        BaseLib::writeValueBinary(fin, v.array()(k));
                    for (size_t k=0; k<n_dummy; k++)
                        BaseLib::writeValueBinary<double>(fin, 0);
                }
            }
        } else {
            offset += bin_data_length + SIZE_OF_BLOCK_LENGTH_TAG;
        }

        if (!_useBinary || !output_data)
            writeDataArrayFooter(fin);
    }

    return true;
}


bool VtuWriter::writeElementGroupID(std::fstream &fin,
                             bool output_data,
                             MeshLib::Mesh& msh,
                             long &offset)
{
    std::string str_format;
    if (!_useBinary)
        str_format = "ascii";
    else
        str_format = "appended";


    if (!_useBinary || !output_data)
        writeDataArrayHeader(fin, this->type_Int, "MatGroup", 0, str_format, offset);
    if (output_data)
    {
        const MeshLib::Element* ele = NULL;
        if (!_useBinary)
        {
            fin << "          ";
            for(long i = 0; i < (long)msh.getNElements(); i++)
            {
                ele = msh.getElement(i);
                fin << ele->getValue() << " ";
            }
            fin << std::endl;
        }
        else
        {
            BaseLib::writeValueBinary<unsigned int>(
                    fin,
                    sizeof(int) *
                    (long)msh.getNElements());
            for (long i = 0; i < (long)msh.getNElements(); i++)
                BaseLib::writeValueBinary(fin, msh.getElement(i)->getValue());
        }
    }
    else
    {
        offset += (long)msh.getNElements() * sizeof(int) +
                  SIZE_OF_BLOCK_LENGTH_TAG;
    }

    if (!_useBinary || !output_data)
        writeDataArrayFooter(fin);


    return true;
}

unsigned char VtuWriter::getVTKCellType(const MeshElemType ele_type)
{
    unsigned char cell_type = 0;

    switch(ele_type)
    {
    case MeshElemType::LINE:               // vtk_line=3
        cell_type = 3;
        break;
    case MeshElemType::QUAD:               // quadrilateral=9
        cell_type = 9;
        break;
    case MeshElemType::HEXAHEDRON:         // hexahedron=12
        cell_type = 12;
        break;
    case MeshElemType::TRIANGLE:           // triangle=5
        cell_type = 5;
        break;
    case MeshElemType::TETRAHEDRON:        // tetrahedron=10
        cell_type = 10;
        break;
    case MeshElemType::PRISM:              // wedge=13
        cell_type = 13;
        break;
    case MeshElemType::PYRAMID:              // pyramid=14
        cell_type = 14;
        break;
    default:
        std::cerr << "***ERROR: NO CORRESPONDING VTK CELL TYPE FOUND. (ELEMENT TYPE=" <<
        MeshLib::MeshElemType2String(ele_type) << ")" << std::endl;
        break;
    }
    return cell_type;
}

void VtuWriter::initialize()
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

size_t VtuWriter::getSizeOfVtkDataType(DataType data_type) const
{
    size_t n = 0;
    switch (data_type) {
    case VtuWriter::Char: n = sizeof(unsigned char); break;
    case VtuWriter::Int: n = sizeof(long); break;
    case VtuWriter::Real: n = sizeof(double); break;
    }

    return n;
}

bool VtuWriter::writeDataArrayHeader(std::fstream &fin,
                                VTK_XML_DATA_TYPE data_type,
                                const std::string &str_name,
                                int nr_components,
                                const std::string &str_format,
                                long offset)
{
    std::string str_data_type;
    switch (data_type)
    {
    case Int8: str_data_type = "Int8";
        break;
    case UInt8: str_data_type = "UInt8";
        break;
    case Int16: str_data_type = "Int16";
        break;
    case UInt16: str_data_type = "UInt16";
        break;
    case Int32: str_data_type = "Int32";
        break;
    case UInt32: str_data_type = "UInt32";
        break;
    case Int64: str_data_type = "Int64";
        break;
    case UInt64: str_data_type = "UInt64";
        break;
    case Float32: str_data_type = "Float32";
        break;
    case Float64: str_data_type = "Float64";
        break;
    }
    fin << "        <DataArray type=\"" << str_data_type << "\"";
    if (str_name != "")
        fin << " Name=\"" << str_name << "\"";
    if (nr_components > 1)
        fin << " NumberOfComponents=\"" << nr_components << "\"";
    fin << " format=\"" << str_format << "\"";
    if (_useBinary)
        fin << " offset=\"" << offset << "\" /";
    fin << ">" << std::endl;

    return true;
}

bool VtuWriter::writeDataArrayFooter(std::fstream &fin)
{
    if (!_useBinary)
        fin << "        </DataArray>" << std::endl;

    return true;
}



bool VtuWriter::writePoints(std::fstream &fin, bool output_data, Mesh& msh, long &offset)
{
    const size_t n_msh_nodes = msh.getNNodes();
    if (output_data)
    {
        if (!_useBinary)
            for (size_t i = 0; i < n_msh_nodes; i++)
            {
                auto* pnt = msh.getNode(i);
                fin << "          " << (*pnt)[0] << " " << (*pnt)[1] << " " << (*pnt)[2] <<
                        std::endl;
            }
        else
        {
            BaseLib::writeValueBinary<unsigned int>(fin, sizeof(double) * 3 * n_msh_nodes);
            for (size_t i = 0; i < n_msh_nodes; i++)
            {
                auto* pnt = msh.getNode(i);
                BaseLib::writeValueBinary(fin, (*pnt)[0]);
                BaseLib::writeValueBinary(fin, (*pnt)[1]);
                BaseLib::writeValueBinary(fin, (*pnt)[2]);
            }
        }
    }
    else if (_useBinary)
        offset += n_msh_nodes * sizeof(double) * 3 + SIZE_OF_BLOCK_LENGTH_TAG;

    return true;
}

bool VtuWriter::writeCellConnectivity(std::fstream &fin,
                                        bool output_data,
                                        Mesh& msh,
                                        long &offset,
                                        long &sum_ele_components)
{
    if (output_data)
    {
        const MeshLib::Element* ele = NULL;
        if (!_useBinary)
            for (long i = 0; i < (long)msh.getNElements(); i++)
            {
                ele = msh.getElement(i);
                fin << "          ";
                for (size_t j = 0; j < ele->getNNodes(); j++)
                    fin << ele->getNodeIndex(j) << " ";
                fin << std::endl;
            }
        else
        {
            BaseLib::writeValueBinary<unsigned int>(fin, sizeof(long) * sum_ele_components);
            for (long i = 0; i < (long)msh.getNElements(); i++)
            {
                ele = msh.getElement(i);
                for (size_t j = 0; j < ele->getNNodes(); j++)
                    BaseLib::writeValueBinary<long>(fin, ele->getNodeIndex(j));
            }
        }
    }
    else if (_useBinary)
    {
        sum_ele_components = 0;
        const MeshLib::Element* ele = NULL;
        for (size_t i = 0; i < msh.getNElements(); i++) {
            ele = msh.getElement(i);
            sum_ele_components += ele->getNNodes();
        }
        offset += sum_ele_components * sizeof(long) + SIZE_OF_BLOCK_LENGTH_TAG;
    }

    return true;
}

bool VtuWriter::writeCellOffset(std::fstream &fin, bool output_data, Mesh& msh, long &offset)
{
    if (output_data)
    {
        const MeshLib::Element* ele = NULL;

        if (!_useBinary)
        {
            fin << "          ";
            long ele_offset = 0;
            for (long i = 0; i < (long)msh.getNElements(); i++)
            {
                ele = msh.getElement(i);
                ele_offset += ele->getNNodes();
                fin << ele_offset << " ";
            }
            fin << std::endl;
        }
        else
        {
            BaseLib::writeValueBinary<unsigned int>(fin,
                                             sizeof(long) * (long)msh.getNElements());
            long ele_offset = 0;
            for (long i = 0; i < (long)msh.getNElements(); i++)
            {
                ele = msh.getElement(i);
                ele_offset += ele->getNNodes();
                BaseLib::writeValueBinary(fin, ele_offset);
            }
        }
    }
    else if (_useBinary)
    {
        offset += (long)msh.getNElements() * sizeof(long) + SIZE_OF_BLOCK_LENGTH_TAG;
    }

    return true;
}

bool VtuWriter::writeCellType(std::fstream &fin, bool output_data, Mesh& msh, long &offset)
{
    if (output_data)
    {
        const MeshLib::Element* ele = NULL;
        if (!_useBinary)
        {
            fin << "          ";
            for(long i = 0; i < (long)msh.getNElements(); i++)
            {
                ele = msh.getElement(i);
                fin << (int)this->getVTKCellType(ele->getGeomType()) << " ";
            }
            fin << std::endl;
        }
        else
        {
            BaseLib::writeValueBinary<unsigned int>(
                    fin,
                    sizeof(unsigned char) *
                    (long)msh.getNElements());
            for(long i = 0; i < (long)msh.getNElements(); i++)
            {
                ele = msh.getElement(i);
                BaseLib::writeValueBinary(fin, this->getVTKCellType(ele->getGeomType()));
            }
        }
    }
    else if (_useBinary)
    {
        offset += (long)msh.getNElements() * sizeof(unsigned char) +
                  SIZE_OF_BLOCK_LENGTH_TAG;
    }

    return true;
}

