/**
 * \author Karsten Rink
 * \date   2010-08-25
 * \brief  Implementation of the project data class.
 *
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ProjectData.h"

#include <algorithm>

// ThirdParty/logog
#include "logog/include/logog.hpp"

#include "BaseLib/FileTools.h"
#include "BaseLib/uniqueInsert.h"

#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"
#include "MeshLib/Mesh.h"

#include "NumLib/ODESolver/TimeDiscretizationBuilder.h"

// FileIO
#include "GeoLib/IO/XmlIO/Boost/BoostXmlGmlInterface.h"
#include "MeshLib/IO/readMeshFromFile.h"

#include "BaseLib/ConfigTree.h"

#include "UncoupledProcessesTimeLoop.h"

#include "ProcessLib/GroundwaterFlow/GroundwaterFlowProcess-fwd.h"
#include "ProcessLib/TES/TESProcess.h"


namespace detail
{
static
void readGeometry(std::string const& fname, GeoLib::GEOObjects & geo_objects)
{
    DBUG("Reading geometry file \'%s\'.", fname.c_str());
    GeoLib::IO::BoostXmlGmlInterface gml_reader(geo_objects);
    gml_reader.readFile(fname);
}
}

ProjectData::ProjectData() = default;

ProjectData::ProjectData(BaseLib::ConfigTree const& project_config,
    std::string const& project_directory,
    std::string const& output_directory)
{
    std::string const geometry_file = BaseLib::copyPathToFileName(
            //! \ogs_file_param{prj__geometry}
            project_config.getConfigParameter<std::string>("geometry"), project_directory
        );
    detail::readGeometry(geometry_file, *_geoObjects);

    std::string const mesh_file = BaseLib::copyPathToFileName(
            //! \ogs_file_param{prj__mesh}
            project_config.getConfigParameter<std::string>("mesh"), project_directory
        );

    MeshLib::Mesh* const mesh = MeshLib::IO::readMeshFromFile(mesh_file);
    if (!mesh) {
        OGS_FATAL("Could not read mesh from \'%s\' file. No mesh added.",
            mesh_file.c_str());
    }
    _mesh_vec.push_back(mesh);

    //! \ogs_file_param{prj__curves}
    parseCurves(project_config.getConfigSubtreeOptional("curves"));

    //! \ogs_file_param{prj__process_variables}
    parseProcessVariables(project_config.getConfigSubtree("process_variables"));

    //! \ogs_file_param{prj__parameters}
    parseParameters(project_config.getConfigSubtree("parameters"));

    //! \ogs_file_param{prj__processes}
    parseProcesses(project_config.getConfigSubtree("processes"));

    //! \ogs_file_param{prj__output}
    parseOutput(project_config.getConfigSubtree("output"), output_directory);

    //! \ogs_file_param{prj__time_stepping}
    parseTimeStepping(project_config.getConfigSubtree("time_stepping"));

    //! \ogs_file_param{prj__linear_solvers}
    parseLinearSolvers(project_config.getConfigSubtree("linear_solvers"));

    //! \ogs_file_param{prj__nonlinear_solvers}
    parseNonlinearSolvers(project_config.getConfigSubtree("nonlinear_solvers"));
}

ProjectData::~ProjectData()
{
    delete _geoObjects;

    for (MeshLib::Mesh* m : _mesh_vec)
        delete m;
}

void ProjectData::addMesh(MeshLib::Mesh* mesh)
{
    std::string name = mesh->getName();
    isMeshNameUniqueAndProvideUniqueName(name);
    mesh->setName(name);
    _mesh_vec.push_back(mesh);
}

std::vector<MeshLib::Mesh*>::const_iterator ProjectData::findMeshByName(
        std::string const& name) const
{
    return const_cast<ProjectData&>(*this).findMeshByName(name);
}

std::vector<MeshLib::Mesh*>::iterator ProjectData::findMeshByName(
        std::string const& name)
{
    return std::find_if(_mesh_vec.begin(), _mesh_vec.end(),
            [&name](MeshLib::Mesh* mesh)
            {
                return mesh && (name == mesh->getName());
            });
}

const MeshLib::Mesh* ProjectData::getMesh(const std::string &name) const
{
    std::vector<MeshLib::Mesh*>::const_iterator it = findMeshByName(name);
    return (it == _mesh_vec.end() ? nullptr : *it);
}

bool ProjectData::removeMesh(const std::string &name)
{
    bool mesh_found = false;
    std::vector<MeshLib::Mesh*>::iterator it = findMeshByName(name);
    while (it != _mesh_vec.end())
    {
        delete *it;
        *it = nullptr;
        it = findMeshByName(name);
        mesh_found = true;
    }

    _mesh_vec.erase(std::remove(_mesh_vec.begin(), _mesh_vec.end(), nullptr),
            _mesh_vec.end());
    return mesh_found;
}

void ProjectData::buildProcesses()
{
    for (auto const& pc : _process_configs)
    {
        //! \ogs_file_param{process__type}
        auto const type = pc.peekConfigParameter<std::string>("type");

        //! \ogs_file_param{process__nonlinear_solver}
        auto const nl_slv_name = pc.getConfigParameter<std::string>("nonlinear_solver");
        auto& nl_slv = BaseLib::getOrError(_nonlinear_solvers, nl_slv_name,
            "A nonlinear solver with the given name has not been defined.");

        auto time_disc = NumLib::createTimeDiscretization<GlobalVector>(
                //! \ogs_file_param{process__time_discretization}
                pc.getConfigSubtree("time_discretization")
            );

        if (type == "GROUNDWATER_FLOW")
        {
            // The existence check of the in the configuration referenced
            // process variables is checked in the physical process.
            // TODO at the moment we have only one mesh, later there can be
            // several meshes. Then we have to assign the referenced mesh
            // here.
            _processes.emplace_back(
                ProcessLib::GroundwaterFlow::
                createGroundwaterFlowProcess<GlobalSetupType>(
                    *_mesh_vec[0], *nl_slv, std::move(time_disc),
                    _process_variables, _parameters, pc));
        }
        else if (type == "TES")
        {
            _processes.emplace_back(
                ProcessLib::TES::createTESProcess<GlobalSetupType>(
                    *_mesh_vec[0], *nl_slv, std::move(time_disc),
                    _process_variables, _parameters, pc));
        }
        else
        {
            OGS_FATAL("Unknown process type: %s", type.c_str());
        }
    }

    // process configs are not needed anymore, so clear the storage
    // in order to trigger config tree checks
    _process_configs.clear();
}

bool ProjectData::meshExists(const std::string &name) const
{
    return findMeshByName(name) != _mesh_vec.end();
}

bool ProjectData::isMeshNameUniqueAndProvideUniqueName(std::string &name) const
{
    int count(0);
    bool isUnique(false);
    std::string cpName;

    while (!isUnique)
    {
        isUnique = true;
        cpName = name;

        count++;
        // If the original name already exists we start to add numbers to name for
        // as long as it takes to make the name unique.
        if (count > 1)
            cpName = cpName + "-" + std::to_string(count);

        for (std::vector<MeshLib::Mesh*>::const_iterator it = _mesh_vec.begin();
                it != _mesh_vec.end(); ++it)
            if ( cpName.compare((*it)->getName()) == 0 )
                isUnique = false;
    }

    // At this point cpName is a unique name and isUnique is true.
    // If cpName is not the original name, "name" is changed and isUnique is set to false,
    // indicating that a vector with the original name already exists.
    if (count > 1)
    {
        isUnique = false;
        name = cpName;
    }
    return isUnique;
}

void ProjectData::parseProcessVariables(
    BaseLib::ConfigTree const& process_variables_config)
{
    DBUG("Parse process variables:")
    if (_geoObjects == nullptr) {
        ERR("Geometric objects are required to define process variables.");
        ERR("No geometric objects present.");
        return;
    }

    // TODO at the moment we have only one mesh, later there
    // can be several meshes. Then we have to check for correct mesh here and
    // assign the referenced mesh below.
    if (_mesh_vec.empty() || _mesh_vec[0] == nullptr) {
        ERR("A mesh is required to define process variables.");
        return;
    }

    // _process_variables.reserve(process_variables_config.size());

    for (auto var_config
         //! \ogs_file_param{prj__process_variables__process_variable}
         : process_variables_config.getConfigSubtreeList("process_variable")) {
        // TODO Extend to referenced meshes.
        _process_variables.emplace_back(var_config, *_mesh_vec[0], *_geoObjects);
    }
}

void ProjectData::parseParameters(BaseLib::ConfigTree const& parameters_config)
{
    using namespace ProcessLib;

    DBUG("Reading parameters:");
    //! \ogs_file_param{prj__parameters__parameter}
    for (auto parameter_config : parameters_config.getConfigSubtreeList("parameter"))
    {
        //! \ogs_file_param{parameter__name}
        auto name = parameter_config.getConfigParameter<std::string>("name");
        //! \ogs_file_param{parameter__type}
        auto type = parameter_config.peekConfigParameter<std::string>("type");

        // Create parameter based on the provided type.
        if (type == "Constant")
        {
            INFO("ConstantParameter: %s.", name.c_str());
            _parameters.push_back(createConstParameter(parameter_config));
            _parameters.back()->name = name;
        }
        else if (type == "MeshProperty")
        {
            INFO("MeshPropertyParameter: %s", name.c_str());
            _parameters.push_back(
                createMeshPropertyParameter(parameter_config, *_mesh_vec[0]));
            _parameters.back()->name = name;
        }
        else
        {
            OGS_FATAL("Cannot construct property of given type \'%s\'.",
                type.c_str());
        }
    }
}

void ProjectData::parseProcesses(BaseLib::ConfigTree const& processes_config)
{
    DBUG("Reading processes:");
    //! \ogs_file_param{prj__processes__process}
    for (auto process_config : processes_config.getConfigSubtreeList("process")) {
        // process type must be specified.
        //! \ogs_file_param{process__type}
        process_config.peekConfigParameter<std::string>("type");
        process_config.ignoreConfigParameter("name");
        _process_configs.push_back(std::move(process_config));
    }
}

void ProjectData::parseOutput(BaseLib::ConfigTree const& output_config,
    std::string const& output_directory)
{
    //! \ogs_file_param{prj__output__type}
    output_config.checkConfigParameter("type", "VTK");
    DBUG("Parse output configuration:");

    _output = ProcessLib::Output<GlobalSetupType>::newInstance(output_config, output_directory);
}

void ProjectData::parseTimeStepping(BaseLib::ConfigTree const& timestepping_config)
{
    DBUG("Reading time loop configuration.");

    _time_loop = ApplicationsLib::createUncoupledProcessesTimeLoop<
        GlobalMatrix, GlobalVector>(timestepping_config);

    if (!_time_loop)
    {
        OGS_FATAL("Initialization of time loop failed.");
    }
}

void ProjectData::parseLinearSolvers(BaseLib::ConfigTree const& config)
{
    DBUG("Reading linear solver configuration.");

    //! \ogs_file_param{prj__linear_solvers__linear_solver}
    for (auto conf : config.getConfigSubtreeList("linear_solver"))
    {
        //! \ogs_file_param{prj__linear_solvers__linear_solver__name}
        auto const name = conf.getConfigParameter<std::string>("name");
        BaseLib::insertIfKeyUniqueElseError(_linear_solvers,
            name,
            MathLib::createLinearSolver<GlobalMatrix, GlobalVector,
                GlobalSetupType::LinearSolver>(&conf),
            "The linear solver name is not unique");
    }
}

void ProjectData::parseNonlinearSolvers(BaseLib::ConfigTree const& config)
{
    DBUG("Reading linear solver configuration.");

    //! \ogs_file_param{prj__nonlinear_solvers__nonlinear_solver}
    for (auto conf : config.getConfigSubtreeList("nonlinear_solver"))
    {
        //! \ogs_file_param{prj__nonlinear_solvers__nonlinear_solver__linear_solver}
        auto const ls_name = conf.getConfigParameter<std::string>("linear_solver");
        auto& linear_solver = BaseLib::getOrError(_linear_solvers,
            ls_name, "A linear solver with the given name does not exist.");

        //! \ogs_file_param{prj__nonlinear_solvers__nonlinear_solver__name}
        auto const name = conf.getConfigParameter<std::string>("name");
        BaseLib::insertIfKeyUniqueElseError(_nonlinear_solvers,
            name,
            NumLib::createNonlinearSolver<GlobalMatrix, GlobalVector>(
                *linear_solver, conf).first,
            "The nonlinear solver name is not unique");
    }
}

static std::unique_ptr<MathLib::PiecewiseLinearInterpolation>
createPiecewiseLinearInterpolation(BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{prj__curves__curve__coords}
    auto coords = config.getConfigParameter<std::vector<double>>("coords");
    //! \ogs_file_param{prj__curves__curve__values}
    auto values = config.getConfigParameter<std::vector<double>>("values");
    if (coords.empty() || values.empty())
    {
        OGS_FATAL("The given co-ordinates or values vector is empty.");
    }
    if (coords.size() != values.size())
    {
        OGS_FATAL("The given co-ordinates and values vector sizes are different.");
    }

    return std::unique_ptr<MathLib::PiecewiseLinearInterpolation>{
        new MathLib::PiecewiseLinearInterpolation{std::move(coords),
                                                  std::move(values)}};
}

void ProjectData::parseCurves(
    boost::optional<BaseLib::ConfigTree> const& config)
{
    if (!config) return;

    DBUG("Reading curves configuration.");

    //! \ogs_file_param{prj__curves__curve}
    for (auto conf : config->getConfigSubtreeList("curve"))
    {
        //! \ogs_file_param{prj__curves__curve__name}
        auto const name = conf.getConfigParameter<std::string>("name");
        BaseLib::insertIfKeyUniqueElseError(
            _curves,
            name,
            createPiecewiseLinearInterpolation(conf),
            "The curve name is not unique.");
    }
}
