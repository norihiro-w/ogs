/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "FunctionLiquidPressure.h"

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

#include "GWTools.h"
#include "FeLiquidFlowAssembler.h"


// local assembler

bool FunctionLiquidPressure::initialize(const boost::property_tree::ptree& option)
{
//	boost::property_tree::write_xml(std::cout, option);
//	std::cout << std::endl;
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
    auto local_assembler = new FeLiquidFlowAssembler(*_feObjects, msh->getCoordinateSystem());

    // set up problem
    _problem = new SolutionLib::FemIVBVProblem(*msh, local_assembler);
    _problem->setTimeSteppingFunction(tim);
    // set up variable
    auto* pressure = _problem->addVariable("pressure"); //internal name
    SolutionLib::FemVariableBuilder varBuilder;
    varBuilder.doit(this->getOutputParameterName(Pressure), option, msh, femData->geo, femData->geo_unique_name, _feObjects, pressure);

    // set up solution
    _solution = new SolutionLib::SingleStepFEM(msh, _problem, femData->linalg_type);
    auto* linear_solver = _solution->getLinearEquationSolver();
    auto optNum = option.get_child("Numerics");
    linear_solver->setOption(optNum);
    _solution->getNonlinearSolver()->setOption(optNum);

    // set initial output
    OutputVariableInfo var(this->getOutputParameterName(Pressure), msh->getID(), OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _solution->getCurrentSolution(0));
    femData->outController.setOutput(var.name, var);

    // initial output parameter
    this->setOutput(Pressure, _solution->getCurrentSolution(pressure->getID()));

    // setup velocity
    _vel = new NumLib::FEMIntegrationPointFunctionVector();
    _vel->initialize(msh);
    _vel_3d = new NumLib::TXWrapped3DVectorFunction<NumLib::FEMIntegrationPointFunctionVector>(_vel, msh->getCoordinateSystem());
    this->setOutput(Velocity, _vel_3d);
    OutputVariableInfo var_v(this->getOutputParameterName(Velocity), msh->getID(), OutputVariableInfo::Element, OutputVariableInfo::Real, 3, _vel);
    femData->outController.setOutput(var_v.name, var_v);

    if (femData->fe_mesh_data[msh_id].node_values.count("PRESSURE1")==0) {
        femData->fe_mesh_data[msh_id].node_values.insert("PRESSURE0", std::vector<double>(msh->getNNodes()));
        femData->fe_mesh_data[msh_id].node_values.insert("PRESSURE1", std::vector<double>(msh->getNNodes()));
    }
    if (femData->fe_mesh_data[msh_id].gp_values.count("VELOCITY1")==0) {
        femData->fe_mesh_data[msh_id].gp_values.insert("VELOCITY0", std::vector<std::vector<MathLib::LocalVector>>(msh->getNElements()));
        femData->fe_mesh_data[msh_id].gp_values.insert("VELOCITY1", std::vector<std::vector<MathLib::LocalVector>>(msh->getNElements()));
    }
    return true;
}

void FunctionLiquidPressure::postSolutionAlgorithm(const NumLib::TimeStep &/*time*/)
{
    INFO("Calculating Darcy velocity within elements from fluid pressure...");
    getDarcyVelocity(*_solution->getCurrentSolution(0), *_vel);
}

void FunctionLiquidPressure::updateOutputParameter(const NumLib::TimeStep &/*time*/)
{
    setOutput(Pressure, _solution->getCurrentSolution(0));
    setOutput(Velocity, _vel_3d);
    auto &msh = _problem->getMesh();
    auto msh_id = msh.getID();
    auto &mshData = THMCLib::Ogs6FemData::getInstance()->fe_mesh_data[msh_id];
    auto &vec_p0 = boost::any_cast<std::vector<double>&>(mshData.node_values.find("PRESSURE0")->second);
    auto &vec_p1 = boost::any_cast<std::vector<double>&>(mshData.node_values.find("PRESSURE1")->second);
    for (std::size_t i=0; i<msh.getNNodes(); i++) {
        vec_p0[i] = vec_p1[i];
        vec_p1[i] = _solution->getCurrentSolution(0)->getValue(i);
    }
    typedef std::vector<std::vector<MathLib::LocalVector>> GpVector;
    auto &vec_v0 = boost::any_cast<GpVector&>(mshData.gp_values.find("VELOCITY0")->second);
    auto &vec_v1 = boost::any_cast<GpVector&>(mshData.gp_values.find("VELOCITY1")->second);
    for (std::size_t i=0; i<msh.getNElements(); i++) {
        for (std::size_t ip=0; ip<_vel->getIntegrationPointValues(i).size(); ip++) {
            vec_v0[i][ip] = vec_v1[i][ip];
            vec_v1[i][ip] = _vel->getIntegrationPointValues(i)[ip];
        }
    }
}

void FunctionLiquidPressure::output(const NumLib::TimeStep &/*time*/)
{
    //update data for output
    auto* femData = THMCLib::Ogs6FemData::getInstance();
    OutputVariableInfo var(this->getOutputParameterName(Pressure), _problem->getMesh().getID(), OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _solution->getCurrentSolution(0));
    femData->outController.setOutput(var.name, var);
    OutputVariableInfo var_v(this->getOutputParameterName(Velocity), _problem->getMesh().getID(), OutputVariableInfo::Element, OutputVariableInfo::Real, 3, _vel);
    femData->outController.setOutput(var_v.name, var_v);
};
