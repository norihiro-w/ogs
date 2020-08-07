/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "BoundaryConditionBuilder.h"

#include "MeshGeoToolsLib/BoundaryElementsSearcher.h"
#include "MeshGeoToolsLib/MeshNodeSearcher.h"
#include "MeshGeoToolsLib/SearchLength.h"
#include "MeshLib/Mesh.h"
#include "ProcessLib/BoundaryCondition/BoundaryConditionConfig.h"
#include "ProcessLib/BoundaryCondition/DirichletBoundaryCondition.h"

#ifdef OGS_USE_PYTHON
#include "Python/PythonBoundaryCondition.h"
#endif

#include "NeumannBoundaryCondition.h"

namespace ProcessLib
{
namespace LIE
{
std::unique_ptr<BoundaryCondition>
BoundaryConditionBuilder::createBoundaryCondition(
    const BoundaryConditionConfig& config,
    const NumLib::LocalToGlobalIndexMap& dof_table,
    const MeshLib::Mesh& bulk_mesh, const int variable_id,
    const unsigned integration_order, const unsigned shapefunction_order,
    const std::vector<std::unique_ptr<ParameterLib::ParameterBase>>& parameters,
    const Process& process)
{
    // Surface mesh and bulk mesh must have equal axial symmetry flags!
    if (config.boundary_mesh.isAxiallySymmetric() !=
        bulk_mesh.isAxiallySymmetric())
    {
        OGS_FATAL(
            "The boundary mesh %s axially symmetric but the bulk mesh %s. Both "
            "must have an equal axial symmetry property.",
            config.boundary_mesh.isAxiallySymmetric() ? "is" : "is not",
            bulk_mesh.isAxiallySymmetric() ? "is" : "is not");
    }

    //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition__type}
    auto const type = config.config.peekConfigParameter<std::string>("type");

    if (type == "Dirichlet")
    {
        return ProcessLib::createDirichletBoundaryCondition(
            config.config, config.boundary_mesh, dof_table, variable_id,
            *config.component_id, parameters);
    }
    if (type == "Neumann")
    {
        return ProcessLib::LIE::createNeumannBoundaryCondition(
            config.config, config.boundary_mesh, dof_table, variable_id,
            *config.component_id, integration_order, shapefunction_order,
            bulk_mesh.getDimension(), parameters);
    }

    if (type == "Python")
    {
#ifdef OGS_USE_PYTHON
        return ProcessLib::createPythonBoundaryCondition(
            config.config, config.boundary_mesh, dof_table, bulk_mesh.getID(),
            variable_id, *config.component_id, integration_order,
            shapefunction_order, bulk_mesh.getDimension());
#else
        OGS_FATAL("OpenGeoSys has not been built with Python support.");
#endif
    }

    OGS_FATAL("Unknown boundary condition type: `%s'.", type.c_str());
}

}  // namespace LIE
}  // namespace ProcessLib
