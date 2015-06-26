/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file VtuWriter.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */


#pragma once

#include <string>
#include <vector>

#include "MeshLib/Mesh.h"
#include "MeshLib/MeshEnums.h"
#include "NumLib/Function/ITXFunction.h"

class VtuWriter
{
public:
    enum DataType {Char, Int, Real};
    enum VTK_XML_DATA_TYPE { Int8, UInt8, Int16, UInt16, Int32, UInt32, Int64, UInt64, Float32,
                         Float64 };

    struct AttributeInfo
    {
        std::string attribute_name;
        DataType data_type;
        size_t nr_of_components;
        NumLib::ITXFunction* f;
        VTK_XML_DATA_TYPE vtk_data_type;

        AttributeInfo(const std::string &name, DataType dtype, size_t n_comp, NumLib::ITXFunction* values)
            : attribute_name(name), data_type(dtype), nr_of_components(n_comp), f(values), vtk_data_type(Float64) 
        {};
        AttributeInfo(const std::string &name, size_t n_comp, NumLib::ITXFunction* values)
            : attribute_name(name), data_type(Real), nr_of_components(n_comp), f(values), vtk_data_type(Float64) 
        {};
    };

    typedef std::pair<std::string, AttributeInfo> PointData;
    typedef std::pair<std::string, AttributeInfo> CellData;

    explicit VtuWriter(bool binary_mode);
    virtual ~VtuWriter(void){}

public:
    bool write(const std::string &vtkfile,
                                  MeshLib::Mesh &msh, std::vector<PointData> &nodal, std::vector<CellData> &elemental, bool outGroupID);

protected:
    void initialize();

    std::string getFormatType() const 
    {
        if (!this->_useBinary) return "ascii";
        else return "appended";
    };

    size_t getSizeOfVtkDataType(DataType data_type) const;

    unsigned char getVTKCellType(const MeshLib::MeshElemType ele_type);

    bool writeDataArrayHeader(std::fstream &fin,
                              VTK_XML_DATA_TYPE data_type,
                              const std::string &str_name,
                              int nr_components,
                              const std::string &str_format,
                              long offset = -1);

    bool writeDataArrayFooter(std::fstream &fin);

    bool writePoints(std::fstream &fin,
                               bool output_data,
                               MeshLib::Mesh& m_msh,
                               long &offset);

    bool writeCellConnectivity(std::fstream &fin,
                                             bool output_data,
                                             MeshLib::Mesh& m_msh,
                                             long &offset,
                                             long &sum_ele_components);

    bool writeCellOffset(std::fstream &fin,
                                       bool output_data,
                                       MeshLib::Mesh& m_msh,
                                       long &offset);

    bool writeCellType(std::fstream &fin,
                                     bool output_data,
                                     MeshLib::Mesh& m_msh,
                                     long &offset);
    
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
    //for binary output
    bool _useBinary;
    bool _isInitialized;
    //Endian(byte order)
    bool _isLittleEndian;
    
    VTK_XML_DATA_TYPE type_UChar;
    VTK_XML_DATA_TYPE type_Int;
    VTK_XML_DATA_TYPE type_UInt;
    VTK_XML_DATA_TYPE type_Long;
    VTK_XML_DATA_TYPE type_Double;
    int SIZE_OF_BLOCK_LENGTH_TAG;
    
};
