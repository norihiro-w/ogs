/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "InitialCondition.h"

#include <boost/optional.hpp>
#include <logog/include/logog.hpp>

#include "BaseLib/ConfigTree.h"
#include "BaseLib/Error.h"

#include "MathLib/Point3d.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/Mesh.h"


namespace ProcessLib
{
std::unique_ptr<InitialCondition> createUniformInitialCondition(
    BaseLib::ConfigTree const& config, int const /*n_components*/)
{
    //! \ogs_file_param{initial_condition__type}
    config.checkConfigParameter("type", "Uniform");

    //! \ogs_file_param{initial_condition__Uniform__value}
    auto value = config.getConfigParameter<double>("value");
    DBUG("Using value %g", value);

    return std::unique_ptr<InitialCondition>(
        new UniformInitialCondition(value));
}

std::unique_ptr<InitialCondition> createMeshPropertyInitialCondition(
    BaseLib::ConfigTree const& config,
    MeshLib::Mesh const& mesh,
    int const n_components)
{
    //! \ogs_file_param{initial_condition__type}
    config.checkConfigParameter("type", "MeshProperty");

    //! \ogs_file_param{initial_condition__MeshProperty__field_name}
    auto field_name = config.getConfigParameter<std::string>("field_name");
    DBUG("Using field_name %s", field_name.c_str());

    if (!mesh.getProperties().hasPropertyVector(field_name))
    {
        OGS_FATAL("The required property %s does not exists in the mesh.",
            field_name.c_str());
    }
    auto const& property =
        mesh.getProperties().template getPropertyVector<double>(field_name);
    if (!property)
    {
        OGS_FATAL("The required property %s is not of the requested type.",
            field_name.c_str());
    }

    if (property->getNumberOfComponents() !=
        static_cast<std::size_t>(n_components))
    {
        OGS_FATAL("The required property %s has different number of components %d, "
            "expected %d.",
            field_name.c_str(), property->getNumberOfComponents(), n_components);
    }
    return std::unique_ptr<InitialCondition>(
        new MeshPropertyInitialCondition(*property));
}

}  // namespace ProcessLib
