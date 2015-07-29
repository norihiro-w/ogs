/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 */

#pragma once

#include <array>
#include <string>
#include <fstream>
#include <vector>

#include "BaseLib/FileTools.h"
#include "GeoLib/Polyline.h"

#include "MeshLib/Mesh.h"
#include "MeshLib/MeshEnums.h"
#include "MeshLib/MeshInformation.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/Node.h"

#include "MeshGeoToolsLib/MeshNodeSearcher.h"
#include "MeshGeoToolsLib/MeshNodesAlongPolyline.h"

#include "NumLib/Function/ITXFunction.h"

class TecplotWriter
{
public:
    TecplotWriter(){}
    virtual ~TecplotWriter(){}

    struct AttributeInfo
    {
        std::string attribute_name;
        size_t nr_of_components;
        NumLib::ITXFunction* f;

        AttributeInfo(const std::string &name, size_t n_comp, NumLib::ITXFunction* values)
            : attribute_name(name), nr_of_components(n_comp), f(values)
        {}
    };

    typedef std::pair<std::string, AttributeInfo> PointData;
    typedef std::pair<std::string, AttributeInfo> CellData;

    void write(const std::string &tec_file_name_base,
               bool init,
               double time_current,
               const MeshLib::Mesh &msh,
               const std::vector<PointData> &nodal,
               const std::vector<CellData> &elemental,
               const GeoLib::GeoObject* geoObject,
               const std::string &geo_obj_name,
               std::size_t geo_id,
               MeshGeoToolsLib::MeshNodeSearcher &nodeSearch);

private:
    std::string getElementType(MeshLib::MeshElemType e_type);

    void writeElementIndeces(const MeshLib::Element &e, std::ostream &os);

    std::string getFileNameSuffix(MeshLib::MeshElemType e_type);

    void writeNodeDataDomain(bool init, double time, const MeshLib::Mesh &msh, const std::vector<PointData> &nodal,
            const std::string &tec_file_name_base);

    void writeElementDataDomain(bool init, double time, const MeshLib::Mesh &msh, const std::vector<CellData> &elemental, const std::string &tec_file_name_base);

    void writeNodeDataPoly(bool init, double time, const MeshLib::Mesh &msh, const std::vector<PointData> &nodal,
            const GeoLib::Polyline &ply, const std::string &ply_name, std::size_t ply_id,
            MeshGeoToolsLib::MeshNodeSearcher &nodeSearch,
            const std::string& tec_file_name_base);

    void writeNodeDataPoint(
            bool init, double time_current, const MeshLib::Mesh &msh, const std::vector<PointData> &nodal,
            const GeoLib::Point &pnt, const std::string &geo_name,
            MeshGeoToolsLib::MeshNodeSearcher &nodeSearch,
            const std::string& tec_file_name_base);

};
