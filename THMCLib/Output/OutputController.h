/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file OutputController.h
 *
 * Created on 2012-07-31 by Norihiro Watanabe
 */

#pragma once

#include <vector>
#include <boost/property_tree/ptree.hpp>

#include "BaseLib/OrderedMap.h"
#include "NumLib/TimeStepping/TimeStep.h"
#include "OutputBuilder.h"
#include "OutputTimingBuilder.h"

namespace GeoLib
{
class GEOObjects;
}

namespace MeshLib
{
class Mesh;
}

namespace ogs6
{

/**
 * \brief Result output controller
 */
class OutputController
{
public:
    ~OutputController()
    {
        clear();
    }

    /**
     *
     * @param option
     * @param output_dir
     * @param project_name
     * @param list_mesh
     * @param geo
     * @param geo_unique_name
     */
    void initialize(boost::property_tree::ptree const& option, const std::string &output_dir,
            const std::string &project_name,
            std::vector<MeshLib::Mesh*> &list_mesh, std::vector<MeshGeoToolsLib::MeshNodeSearcher*> &list_nodeSearcher,
            GeoLib::GEOObjects &geo, const std::string &geo_unique_name);

    /**
     *
     * @param time
     * @return
     */
    bool isActive(const NumLib::TimeStep &time) const;

    /**
     *
     * @param time
     */
    void outputData(const NumLib::TimeStep &time);

    /**
     *
     * @param var_name
     * @param var
     */
    void setOutput(const std::string &var_name, OutputVariableInfo &var)
    {
        _map_name_var.insert(var_name, var);
    }

    void clear()
    {
        //BaseLib::releaseObjectsInStdVector(_list_output);
        _map_name_var.clear();
    }

private:
    std::vector<IOutput*> _list_output;
    BaseLib::OrderedMap<std::string, OutputVariableInfo> _map_name_var;
};

}
