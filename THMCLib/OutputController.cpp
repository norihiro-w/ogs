/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 */

#include "OutputController.h"

#include "GeoLib/GEOObjects.h"
#include "MeshLib/Mesh.h"

namespace ogs6
{

void OutputController::initialize(
        boost::property_tree::ptree const& option, const std::string &output_dir, const std::string &project_name,
        std::vector<MeshLib::Mesh*> &list_mesh, std::vector<MeshGeoToolsLib::MeshNodeSearcher*> &list_nodeSearcher, GeoLib::GEOObjects &geo, const std::string &geo_unique_name)
{
    //Ogs6FemData* femData = Ogs6FemData::getInstance();
    OutputBuilder outBuilder;
    OutputTimingBuilder outTimBuilder;

    auto &opOutList = option.get_child("outputList");
    auto range = opOutList.equal_range("output");
    for (auto it=range.first; it!=range.second; ++it)
    {
        auto& op = it->second;
        std::string data_type = op.get<std::string>("dataType");
        IOutput* out = outBuilder.create(data_type);
        out->setOutputPath(output_dir, project_name);

        size_t msh_id = op.get<size_t>("meshID");
        out->setMesh(list_mesh[msh_id]);
        out->setMeshNodeSearcher(list_nodeSearcher[msh_id]);

        std::string time_type = op.get<std::string>("timeType");
        size_t out_steps = op.get<size_t>("timeSteps");
        out->setOutputTiming(outTimBuilder.create(time_type, out_steps));

        std::vector<std::string> nod_var_name;
        auto range = op.equal_range("nodeValue");
        for (auto it=range.first; it!=range.second; ++it)
        {
            nod_var_name.push_back(it->second.get<std::string>("name"));
        }
        out->addNodalVariable(nod_var_name);
        std::vector<std::string> ele_var_name;
        auto rangeE = op.equal_range("elementValue");
        for (auto it=rangeE.first; it!=rangeE.second; ++it)
        {
            ele_var_name.push_back(it->second.get<std::string>("name"));
        }
        out->addElementalVariable(ele_var_name);

        std::string geo_type_name = op.get<std::string>("geoType");
        std::string geo_name = op.get<std::string>("geoName");
        const GeoLib::GeoObject* geo_obj = geo.getGeoObject(geo_unique_name, GeoLib::convertGeoType(geo_type_name), geo_name);
        auto geotype = GeoLib::convertGeoType(geo_type_name);
        if (geotype == GeoLib::GEOTYPE::POLYLINE) {
            std::size_t ply_id = 0;
            geo.getPolylineVecObj(geo_unique_name)->getElementIDByName(geo_name, ply_id);
            out->setGeometryIndex(ply_id);
        }
        out->setGeometry(geo_obj);
        out->setGeometryName(geo_name);

        _list_output.push_back(out);
    }
}

bool OutputController::isActive(const NumLib::TimeStep &time) const
{
    bool doOutput = false;
    for (size_t i=0; i<_list_output.size(); i++) {
        if (_list_output[i]->isActive(time)) {
            doOutput = true;
            break;
        }
    }
    return doOutput;
}

void OutputController::outputData(const NumLib::TimeStep &time)
{
    if (!isActive(time)) return;

    for (size_t i=0; i<_list_output.size(); i++) {
        if (_list_output[i]->isActive(time)) {
            _list_output[i]->write(time, _map_name_var);
        }
    }
}

}
