/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 */


#pragma once

#include <string>
#include <vector>

#include "MeshLib/Mesh.h"
#include "MeshLib/MeshEnums.h"
#include "NumLib/Function/ITXFunction.h"
#include "VtuWriter.h"

class PVtuWriter
{
public:
//    enum DataType {Char, Int, Real};
//    enum VTK_XML_DATA_TYPE { Int8, UInt8, Int16, UInt16, Int32, UInt32, Int64, UInt64, Float32,
//                         Float64 };

    typedef std::pair<std::string, VtuWriter::AttributeInfo> PointData;
    typedef std::pair<std::string, VtuWriter::AttributeInfo> CellData;

    explicit PVtuWriter();
    virtual ~PVtuWriter(void){}

public:
    bool write(const std::string &vtkfile,
                                  MeshLib::Mesh &msh, std::vector<PointData> &nodal, std::vector<CellData> &elemental, bool outGroupID);

protected:
    void initialize();

    std::string getFormatType() const { return "ascii"; }

    size_t getSizeOfVtkDataType(VtuWriter::DataType data_type) const;

    unsigned char getVTKCellType(const MeshLib::MeshElemType ele_type);

    bool writeDataArrayHeader(std::fstream &fin,
                              VtuWriter::VTK_XML_DATA_TYPE data_type,
                              const std::string &str_name,
                              int nr_components,
                              const std::string &str_format,
                              long offset = -1);

    bool writePointData(std::fstream &fin,
                                bool output_data,
                                std::vector<PointData> &nodal,
                                MeshLib::Mesh &m_msh,
                                long &offset);
    bool writeCellData(std::fstream &fin,
                                  bool output_data,
                                  std::vector<CellData> &elemental,
                                    MeshLib::Mesh &m_msh,
                                  long &offset);

    bool writeElementGroupID(std::fstream &fin,
                                 bool output_data,
                                 MeshLib::Mesh& msh,
                                 long &offset);

protected:
    bool _isInitialized;
    //Endian(byte order)
    bool _isLittleEndian;
    
    VtuWriter::VTK_XML_DATA_TYPE type_UChar;
    VtuWriter::VTK_XML_DATA_TYPE type_Int;
    VtuWriter::VTK_XML_DATA_TYPE type_UInt;
    VtuWriter::VTK_XML_DATA_TYPE type_Long;
    VtuWriter::VTK_XML_DATA_TYPE type_Double;
    int SIZE_OF_BLOCK_LENGTH_TAG;
    
};
