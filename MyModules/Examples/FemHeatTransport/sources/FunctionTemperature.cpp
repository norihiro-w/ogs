/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 */

#include "FunctionTemperature.h"

#include <logog/include/logog.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include "BaseLib/CodingTools.h"

#include "MathLib/LinAlg/ILinearSolver.h"

#include "MeshLib/ElementCoordinatesMappingLocal.h"

#include "NumLib/Function/FemNodalFunction.h"
#include "NumLib/Function/TXFunctionBuilder.h"
#include "NumLib/Fem/Integration/IIntegration.h"
#include "NumLib/Fem/Integration/GaussIntegrationPolicy.h"

#include "FileIO/OutputIO/OutputBuilder.h"
#include "FileIO/OutputIO/OutputTimingBuilder.h"

#include "SolutionLib/Fem/FemSourceTerm.h"

#include "Ogs6FemData.h"
#include "FemVariableBuilder.h"

#include "FeHeatTransportAssembler.h"

bool FunctionTemperature::initialize(const boost::property_tree::ptree &option)
{
    auto* femData = THMCLib::Ogs6FemData::getInstance();
    auto opt_msh_id = option.get_optional<size_t>("MeshID");
    if (!opt_msh_id) {
        ERR("Mesh ID not given for LIQUID_FLOW");
        return false;
    }
    size_t msh_id = *opt_msh_id;
    auto time_id = option.get_optional<size_t>("TimeGroupID");
    NumLib::ITimeStepAlgorithm* tim = nullptr;
    if (time_id)
        tim = femData->list_tim[*time_id];

    //mesh and FE objects
    MeshLib::Mesh* msh = femData->list_mesh[msh_id];
    _feObjects = new NumLib::LagrangeFeObjectContainer(msh);

    // local assemblers
    auto local_assembler = new FeHeatTransportAssembler(*_feObjects, msh->getCoordinateSystem());

    // set up problem
    _problem = new SolutionLib::FemIVBVProblem(*msh, local_assembler);
    _problem->setTimeSteppingFunction(tim);
    // set up variable
    auto* pressure = _problem->addVariable("temperature"); //internal name
    SolutionLib::FemVariableBuilder varBuilder;
    varBuilder.doit(this->getOutputParameterName(Temperature), option, femData->list_nodeSearcher[msh_id], femData->list_beSearcher[msh_id],
                    femData->geo, femData->geo_unique_name, _feObjects, pressure);

    // set up solution
    _solution = new SolutionLib::SingleStepFEM(msh, _problem, femData->linalg_type);
    auto* linear_solver = _solution->getLinearEquationSolver();
    auto optNum = option.get_child("Numerics");
    linear_solver->setOption(optNum);
    _solution->getNonlinearSolver()->setOption(optNum);

    // set initial output
    OutputVariableInfo var(this->getOutputParameterName(Temperature), msh->getID(), OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _solution->getCurrentSolution(0));
    femData->outController.setOutput(var.name, var);

    // initial output parameter
    this->setOutput(Temperature, _solution->getCurrentSolution(pressure->getID()));

    if (femData->fe_mesh_data[msh_id].node_values.count("TEMPERATURE1")==0) {
        femData->fe_mesh_data[msh_id].node_values.insert("TEMPERATURE0", std::vector<double>(msh->getNNodes()));
        femData->fe_mesh_data[msh_id].node_values.insert("TEMPERATURE1", std::vector<double>(msh->getNNodes()));
    }
    return true;
}

void FunctionTemperature::initializeTimeStep(const NumLib::TimeStep &/*time*/)
{
//    const NumLib::ITXFunction *vel = this->getInput<NumLib::ITXFunction>(Velocity);
}

void FunctionTemperature::updateOutputParameter(const NumLib::TimeStep &/*time*/)
{
    setOutput(Temperature, _solution->getCurrentSolution(0));
    auto &msh = _problem->getMesh();
    auto msh_id = msh.getID();
    auto &mshData = THMCLib::Ogs6FemData::getInstance()->fe_mesh_data[msh_id];
    auto &vec_T0 = boost::any_cast<std::vector<double>&>(mshData.node_values.find("TEMPERATURE0")->second);
    auto &vec_T1 = boost::any_cast<std::vector<double>&>(mshData.node_values.find("TEMPERATURE1")->second);
    for (std::size_t i=0; i<msh.getNNodes(); i++) {
        vec_T0[i] = vec_T1[i];
        vec_T1[i] = _solution->getCurrentSolution(0)->getValue(i);
    }
}

void FunctionTemperature::output(const NumLib::TimeStep &/*time*/)
{
    //update data for output
    auto* femData = THMCLib::Ogs6FemData::getInstance();
    OutputVariableInfo var(this->getOutputParameterName(Temperature), _problem->getMesh().getID(), OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _solution->getCurrentSolution(0));
    femData->outController.setOutput(var.name, var);
};
