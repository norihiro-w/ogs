/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ShapeFunctionParameter.h"
#include "BaseLib/ConfigTree.h"
#include "MeshLib/Mesh.h"

namespace ProcessLib
{
std::unique_ptr<ParameterBase> createShapeFunctionParameter(
    std::string const& name, BaseLib::ConfigTree const& config,
    MeshLib::Mesh const& mesh)
{
    //! \ogs_file_param{prj__parameters__parameter__type}
    config.checkConfigParameter("type", "ShapeFunction");
    //! \ogs_file_param{prj__parameters__parameter__MeshNode__field_name}
    auto const field_name = config.getConfigParameter<std::string>("field_name");
    DBUG("Using field_name %s", field_name.c_str());
    auto const shape_function_order = config.getConfigParameter<unsigned>("order");

    // TODO other data types than only double
    auto const& property =
        mesh.getProperties().getPropertyVector<double>(field_name);

    if (property->getMeshItemType() != MeshLib::MeshItemType::Node) {
        OGS_FATAL("The mesh property `%s' is not a nodal property.",
                  field_name.c_str());
    }

    return std::make_unique<ShapeFunctionParameter<double>>(
        name, mesh, shape_function_order, *property);
}

}  // ProcessLib
