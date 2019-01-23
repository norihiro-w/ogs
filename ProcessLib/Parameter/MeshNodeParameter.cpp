/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "MeshNodeParameter.h"
#include "BaseLib/ConfigTree.h"
#include "MeshLib/Mesh.h"

namespace ProcessLib
{
std::unique_ptr<ParameterBase> createMeshNodeParameter(
    std::string const& name, BaseLib::ConfigTree const& config,
    MeshLib::Mesh const& mesh)
{
    //! \ogs_file_param{prj__parameters__parameter__type}
    config.checkConfigParameter("type", "MeshNode");
    //! \ogs_file_param{prj__parameters__parameter__MeshNode__field_name}
    auto const field_name = config.getConfigParameter<std::string>("field_name");
    DBUG("Using field_name %s", field_name.c_str());

    auto const hasProperty =
        mesh.getProperties().existsPropertyVector<double>(field_name);
    auto const create =
    //! \ogs_file_param{prj__parameters__parameter__MeshElement__create}
        config.getConfigParameterOptional<std::string>("create");
    if (!hasProperty && create && create.get() == "true")
    {
        DBUG("Property %s does not exit. create", field_name.c_str());
        auto opt_n_comp =
            //! \ogs_file_param{prj__parameters__parameter__MeshElement__components}
            config.getConfigParameterOptional<unsigned>("components");
        unsigned const n_comp = opt_n_comp ? opt_n_comp.get() : 1;

        auto new_property = const_cast<MeshLib::Mesh&>(mesh)
                                .getProperties()
                                .createNewPropertyVector<double>(
                                    field_name, MeshLib::MeshItemType::Node, n_comp);
        new_property->resize(mesh.getNumberOfNodes() * n_comp);
    }

    // TODO other data types than only double
    auto const& property =
        mesh.getProperties().getPropertyVector<double>(field_name);

    if (property->getMeshItemType() != MeshLib::MeshItemType::Node) {
        OGS_FATAL("The mesh property `%s' is not a nodal property.",
                  field_name.c_str());
    }

    return std::make_unique<MeshNodeParameter<double>>(name, *property);
}

}  // ProcessLib
