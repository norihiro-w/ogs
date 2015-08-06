/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 */

#include "TecplotWriter.h"


void TecplotWriter::write(const std::string &tec_file_name_base,
           bool init,
           double time_current,
           const MeshLib::Mesh &msh,
           const std::vector<PointData> &nodal,
           const std::vector<CellData> &elemental,
           const GeoLib::GeoObject* geoObject,
           const std::string &geo_obj_name,
           std::size_t geo_id,
           MeshGeoToolsLib::MeshNodeSearcher &nodeSearch)
{
    if (geoObject->getGeoType() == GeoLib::GEOTYPE::GEODOMAIN) {
        writeNodeDataDomain(init, time_current, msh, nodal, tec_file_name_base);
        writeElementDataDomain(init, time_current, msh, elemental, tec_file_name_base);
    } else if (geoObject->getGeoType() == GeoLib::GEOTYPE::SURFACE) {
        //TODO
    } else if (geoObject->getGeoType() == GeoLib::GEOTYPE::POLYLINE) {
        GeoLib::Polyline const* const ply (dynamic_cast<GeoLib::Polyline const* const>(geoObject));
        writeNodeDataPoly(init, time_current, msh, nodal, *ply, geo_obj_name, geo_id, nodeSearch, tec_file_name_base);
    } else if (geoObject->getGeoType() == GeoLib::GEOTYPE::POINT) {
        GeoLib::Point const* const pnt (dynamic_cast<GeoLib::Point const* const>(geoObject));
        writeNodeDataPoint(init, time_current, msh, nodal, *pnt, geo_obj_name, nodeSearch, tec_file_name_base);
    }

}

void TecplotWriter::writeNodeDataDomain(bool init, double time, const MeshLib::Mesh &msh, const std::vector<PointData> &nodal,
        const std::string &tec_file_name_base)
{
    if(nodal.size() == 0)
        return;

    const std::array<unsigned, 7> vec_ele_type_n = MeshLib::MeshInformation::getNumberOfElementTypes(msh);
    std::vector<MeshLib::MeshElemType> ele_types;
    ele_types.push_back(MeshLib::MeshElemType::LINE);
    ele_types.push_back(MeshLib::MeshElemType::TRIANGLE);
    ele_types.push_back(MeshLib::MeshElemType::QUAD);
    ele_types.push_back(MeshLib::MeshElemType::TETRAHEDRON);
    ele_types.push_back(MeshLib::MeshElemType::HEXAHEDRON);
    ele_types.push_back(MeshLib::MeshElemType::PYRAMID);
    ele_types.push_back(MeshLib::MeshElemType::PRISM);

    // Output files for each mesh type
    for (std::size_t i = 0; i < vec_ele_type_n.size(); i++)
    {
        if (vec_ele_type_n[i] == 0) continue;

        auto eletype = ele_types[i];
        std::string tec_file_name = tec_file_name_base + "_domain_" + getFileNameSuffix(eletype) + ".tec";
        if (init)
            BaseLib::truncateFile(tec_file_name);

        std::ios_base::openmode omode = std::ios::app;
        if (init) omode = std::ios::trunc;
        std::ofstream tec_file (tec_file_name.data(), omode);
        if (!tec_file.good())
            return;
        tec_file.setf(std::ios::scientific, std::ios::floatfield);
        tec_file.precision(12);

        // header
        tec_file << "VARIABLES  = \"X\",\"Y\",\"Z\"";
        for (auto v : nodal)
            tec_file << ", \"" << v.first << "\"";
        tec_file << "\n";
        tec_file << "ZONE T=\"";
        tec_file << time << "s\", ";
        tec_file << "N=" << msh.getNNodes() << ", ";
        tec_file << "E=" << vec_ele_type_n[i] << ", ";
        tec_file << "F=" << "FEPOINT" << ", ";
        tec_file << "ET=" << getElementType(eletype);
        tec_file << "\n";
        tec_file << "STRANDID=1, SOLUTIONTIME=" << time << "\n";

        // node data
        for (size_t j = 0; j < msh.getNNodes(); j++)
        {
            auto node = msh.getNode(j);

            // XYZ
            for (size_t i = 0; i < 3; i++)
                tec_file <<(*node)[i] << " ";

            for (size_t k = 0; k < nodal.size(); k++)
            {
                MathLib::LocalMatrix v;
                nodal[k].second.f->eval(NumLib::TXPosition(NumLib::TXPosition::Node, j, msh.getNode(j)->getCoords()), v);
#ifdef OGS_USE_EIGEN
                const size_t n_dummy = nodal[k].second.nr_of_components - static_cast<size_t>(v.array().size());
                for (int k=0; k<v.array().size(); k++)
                    tec_file << v.array()(k) << " ";
                for (size_t k=0; k<n_dummy; k++)
                    tec_file << 0.0 << " ";
#endif
            }
            tec_file << "\n";
        }

        // element data
        for (auto e : msh.getElements())
            if (e->getGeomType() == eletype)
                writeElementIndeces(*e, tec_file);

        tec_file.close();
    }
}

void TecplotWriter::writeElementDataDomain(bool init, double time, const MeshLib::Mesh &msh, const std::vector<CellData> &elemental, const std::string &tec_file_name_base)
{
    std::string tec_file_name = tec_file_name_base + "_domain_ele.tec";
    if (init)
        BaseLib::truncateFile(tec_file_name);
    std::ofstream tec_file (tec_file_name.data(), std::ios::app);
    tec_file.setf(std::ios::scientific, std::ios::floatfield);
    tec_file.precision(12);
    if (!tec_file.good())
        return;

    // Write Header I: variables
    tec_file << "VARIABLES = \"X\",\"Y\",\"Z\",\"VX\",\"VY\",\"VZ\"";
    for (size_t i = 0; i < elemental.size(); i++)
        if (elemental[i].first.find("VELOCITY") == std::string::npos)
            tec_file << "," << elemental[i].first;
    tec_file << "\n";

    // Write Header II: zone
    tec_file << "ZONE T=\"";
    tec_file << time << "s\", ";
    tec_file << "I=" << msh.getNElements() << ", ";
    tec_file << "F=POINT" << ", ";
    tec_file << "C=BLACK";
    tec_file << "\n";

    for (size_t i = 0; i < msh.getNElements(); i++)
    {
        auto e = msh.getElement(i);
        auto xyz(e->getCenterOfGravity());
        tec_file << xyz[0] << " " << xyz[1] << " " << xyz[2] << " ";
        for (size_t j = 0; j < elemental.size(); j++)
        {
            MathLib::LocalMatrix v;
            elemental[j].second.f->eval(NumLib::TXPosition(NumLib::TXPosition::Element, i), v);
#ifdef OGS_USE_EIGEN
            const size_t n_dummy = elemental[j].second.nr_of_components - static_cast<size_t>(v.array().size());
            for (int k = 0; k < v.array().size(); k++)
                tec_file << v.array()(k) << " ";
            for (size_t k = 0; k < n_dummy; k++)
                tec_file << 0.0 << " ";
#endif
        }
        tec_file << "\n";
    }
    tec_file << "\n";

    tec_file.close();
}

void TecplotWriter::writeNodeDataPoly(bool init, double time, const MeshLib::Mesh &msh, const std::vector<PointData> &nodal,
        const GeoLib::Polyline &ply, const std::string &ply_name, std::size_t ply_id,
        MeshGeoToolsLib::MeshNodeSearcher &nodeSearch,
        const std::string& tec_file_name_base)
{
    const std::string tec_file_name = tec_file_name_base + "_ply_" + ply_name + + "_t" + std::to_string(ply_id) + ".tec";

    if (init)
        BaseLib::truncateFile(tec_file_name);
    std::ofstream tec_file(tec_file_name.data(), std::ios::app);
    if (!tec_file.good())
        return;

    tec_file.setf(std::ios::scientific, std::ios::floatfield);
    tec_file.precision(12);

    //--------------------------------------------------------------------
    // Write header
//        if (!fileExist)
//        {
        //project_title;
        std::string project_title_string = "Profiles along polylines";
        tec_file << " TITLE = \"" << project_title_string << "\"" << "\n";
        tec_file << " VARIABLES = \"DIST\" ";
        for (size_t k = 0; k < nodal.size(); k++)
        {
            tec_file << "\"" << nodal[k].first << "\" ";
        }
        tec_file << "\n";
//        }
    tec_file << " ZONE T=\"TIME=" << time << "\"" << "\n";

    // Write node values
    auto searchply = nodeSearch.getMeshNodesAlongPolyline(ply);
    auto nodeids = searchply.getNodeIDs();
    auto distances = searchply.getDistOfProjNodeFromPlyStart();
    for (std::size_t i=0; i<nodeids.size(); i++)
    {
        auto nodeid = nodeids[i];
        auto dist = distances[i];
        tec_file << dist << " ";
        for (size_t k = 0; k < nodal.size(); k++)
        {
            MathLib::LocalMatrix v;
            nodal[k].second.f->eval(NumLib::TXPosition(NumLib::TXPosition::Node, nodeid, msh.getNode(nodeid)->getCoords()), v);
#ifdef OGS_USE_EIGEN
            const size_t n_dummy = nodal[k].second.nr_of_components - static_cast<size_t>(v.array().size());
            for (int k=0; k<v.array().size(); k++)
                tec_file << v.array()(k) << " ";
            for (size_t k=0; k<n_dummy; k++)
                tec_file << 0.0 << " ";
#endif
        }
        tec_file << "\n";
    }

    tec_file.close();
}

void TecplotWriter::writeNodeDataPoint(
        bool init, double time_current, const MeshLib::Mesh &msh, const std::vector<PointData> &nodal,
        const GeoLib::Point &pnt, const std::string &geo_name,
        MeshGeoToolsLib::MeshNodeSearcher &nodeSearch,
        const std::string& tec_file_name_base)
{
    const std::string tec_file_name = tec_file_name_base + "_" + geo_name + ".tec";
    if (init)
        BaseLib::truncateFile(tec_file_name);
    std::ofstream tec_file(tec_file_name.data(), std::ios::app);
    if (!tec_file.good())
        return;

    tec_file.setf(std::ios::scientific, std::ios::floatfield);
    tec_file.precision(12);

    const std::string sep = " ";

    // Write header
    if (init)
    {
        const std::string project_title_string ("Time curves in points");
        tec_file << " TITLE = \"" << project_title_string << "\"" << "\n";
        tec_file << " VARIABLES = ";
        tec_file << "\"TIME\"";
        for (size_t k = 0; k < nodal.size(); k++)
            tec_file << sep << "\"" << nodal[k].first << "\"";
        tec_file << "\n";

        tec_file << " ZONE T=\"POINT=" << geo_name << "\"" << "\n";
    }

    // Write data
    tec_file << time_current;

    auto nodeids = nodeSearch.getMeshNodeIDForPoint(pnt);
    for (std::size_t i=0; i<nodeids.size(); i++)
    {
        auto nodeid = nodeids[i];
        for (size_t k = 0; k < nodal.size(); k++)
        {
            MathLib::LocalMatrix v;
            nodal[k].second.f->eval(NumLib::TXPosition(NumLib::TXPosition::Node, nodeid, msh.getNode(nodeid)->getCoords()), v);
#ifdef OGS_USE_EIGEN
            const size_t n_dummy = nodal[k].second.nr_of_components - static_cast<size_t>(v.array().size());
            for (int k=0; k<v.array().size(); k++)
                tec_file << v.array()(k) << " ";
#else
            const size_t n_dummy = nodal[k].second.nr_of_components - static_cast<size_t>(v.rows()*v.columns());
            for (int k=0; k<v.rows(); k++)
                for (int l=0; l<v.columns(); l++)
                    tec_file << v(k,l) << " ";
#endif
            for (size_t k=0; k<n_dummy; k++)
                tec_file << 0.0 << " ";
        }
        tec_file << "\n";
    }
    tec_file << "\n";
    tec_file.close();
}

std::string TecplotWriter::getFileNameSuffix(MeshLib::MeshElemType e_type)
{
    switch (e_type)
    {
    case MeshLib::MeshElemType::LINE:
        return "line";
    case MeshLib::MeshElemType::QUAD:
        return "quad";
    case MeshLib::MeshElemType::HEXAHEDRON:
        return "hex";
    case MeshLib::MeshElemType::TRIANGLE:
        return "tri";
    case MeshLib::MeshElemType::TETRAHEDRON:
        return "tet";
    case MeshLib::MeshElemType::PRISM:
        return "pris";
    case MeshLib::MeshElemType::PYRAMID:
        return "pyra";
    default:
        return "";
    }
}

std::string TecplotWriter::getElementType(MeshLib::MeshElemType e_type)
{
    std::string eleType;
    switch (e_type)
    {
    case MeshLib::MeshElemType::LINE:
    case MeshLib::MeshElemType::QUAD:
    case MeshLib::MeshElemType::TRIANGLE:
        eleType = "QUADRILATERAL";
        break;
    case MeshLib::MeshElemType::TETRAHEDRON:
        eleType = "TETRAHEDRON";
        break;
    case MeshLib::MeshElemType::HEXAHEDRON:
    case MeshLib::MeshElemType::PRISM:
    case MeshLib::MeshElemType::PYRAMID:
        eleType = "BRICK";
        break;
    default:
        eleType = "";
        break;
    }
    return eleType;
}

void TecplotWriter::writeElementIndeces(const MeshLib::Element &e, std::ostream &os)
{
    std::string deli = "  ";
    std::vector<std::size_t> nodes_index;
    for (unsigned i=0; i<e.getNNodes(); i++)
        nodes_index.push_back(e.getNodeIndex(i));
    if(e.getGeomType() == MeshLib::MeshElemType::LINE)
        os << nodes_index[0] + 1 << deli << nodes_index[1] + 1 << deli << nodes_index[1] +
        1 << deli << nodes_index[0] + 1;

    else if(e.getGeomType() == MeshLib::MeshElemType::TRIANGLE)
        os << nodes_index[0] + 1 << deli << nodes_index[1] + 1 << deli << nodes_index[2] +
        1 << deli << nodes_index[0] + 1;
    else if(e.getGeomType() == MeshLib::MeshElemType::PRISM)
        os << nodes_index[0] + 1 << deli << nodes_index[0] + 1 << deli << nodes_index[1] +
        1 << deli << nodes_index[2] + 1 << deli
           << nodes_index[3] + 1 << deli << nodes_index[3] + 1 << deli << nodes_index[4] +
        1 << deli << nodes_index[5] + 1 << deli;
    else if(e.getGeomType() == MeshLib::MeshElemType::PYRAMID)
        os << nodes_index[0] + 1 << deli << nodes_index[1] + 1 << deli << nodes_index[2] +
        1 << deli << nodes_index[3] + 1 << deli
           << nodes_index[4] + 1 << deli << nodes_index[4] + 1 << deli << nodes_index[4] +
        1 << deli << nodes_index[4] + 1 << deli;
    else
        for(int i = 0; i < e.getNNodes(); i++)
            os << nodes_index[i] + 1 << deli;
    os << '\n';
}
