/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 */

#include "TecplotOutput.h"

#include <vector>
#include <map>

#include "BaseLib/CodingTools.h"
#include "NumLib/Function/ITXFunction.h"

#include "FileIO/OutputIO/TecplotWriter.h"

void TecplotOutput::write(const NumLib::TimeStep &current_time,
            BaseLib::OrderedMap<std::string, OutputVariableInfo> &data)
{
    // prepare vtu data
    std::vector<TecplotWriter::PointData> node_values;
    std::vector<TecplotWriter::CellData> ele_values;

    for (BaseLib::OrderedMap<std::string, OutputVariableInfo>::iterator itr = data.begin(); itr!=data.end(); ++itr) {
        // pick up variables for this mesh and specified by users
        if (itr->second.mesh_id != getMesh()->getID()) continue;
        if ((itr->second.object_type == OutputVariableInfo::Node && hasNodalVariable(itr->first))
            || (itr->second.object_type == OutputVariableInfo::Element && hasElementalVariable(itr->first))) {
            OutputVariableInfo &var = itr->second;
            TecplotWriter::AttributeInfo attr(var.name, var.nr_of_components, var.value);
            if (var.object_type==OutputVariableInfo::Node) {
                node_values.push_back(TecplotWriter::PointData(var.name, attr));
            } else if (var.object_type==OutputVariableInfo::Element) {
                ele_values.push_back(TecplotWriter::CellData(var.name, attr));
            }
        }
    }

    if (node_values.size() == 0 && ele_values.size() == 0) {
        WARN("***Warning: Asked to write results but specified data not found. Skip this output. Please check consistency of variable names in process and output settings.");
        return;
    }

    // set base file name for output
    const std::string tec_file_basename = getOutputBaseName();
    std::string tec_file_path = getOutputDir();
    if (tec_file_path.length()>0) tec_file_path += "/";

    // write VTU file
    std::string tec_file_name_base = tec_file_basename;
    std::string tec_file_name = tec_file_path + tec_file_name_base;
    INFOa("Writing results...: %s", tec_file_name.c_str());
    TecplotWriter writer;
    writer.write(tec_file_name_base,
            _init, current_time.current(), *getMesh(), node_values, ele_values,
            getGeometry(), getGeometryName(), getGeometryIndex(), *getMeshNodeSearcher());
    _init = false;
}
