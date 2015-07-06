/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 */

#include "FunctionPressureToElementVelocity.h"

#include "BaseLib/CodingTools.h"
#include "NumLib/Function/FemNodalFunction.h"
#include "THMCLib/Ogs6FemData.h"
#include "GWTools.h"

bool FunctionPressureToElementVelocity::initialize(const boost::property_tree::ptree &option)
{
    auto* femData = THMCLib::Ogs6FemData::getInstance();

    size_t msh_id = option.get<size_t>("MeshID");
    _msh = femData->list_mesh[msh_id];

    _vel = new NumLib::FEMIntegrationPointFunctionVector();
    _vel->initialize(_msh);

    // set initial output
    OutputVariableInfo var(this->getOutputParameterName(Velocity), msh_id, OutputVariableInfo::Element, OutputVariableInfo::Real, 3, _vel);
    femData->outController.setOutput(var.name, var);

    // initial output parameter
    _vel_3d = new NumLib::TXWrapped3DVectorFunction<NumLib::FEMIntegrationPointFunctionVector>(_vel, _msh->getCoordinateSystem());
    this->setOutput(Velocity, _vel_3d);

    return true;
}

void FunctionPressureToElementVelocity::accept(const NumLib::TimeStep &/*time*/)
{
    //std::cout << "Velocity=" << std::endl;
    //_vel->printout();
    //update data for output
    auto* femData = THMCLib::Ogs6FemData::getInstance();
    OutputVariableInfo var(this->getOutputParameterName(Velocity), _msh->getID(), OutputVariableInfo::Element, OutputVariableInfo::Real, 3, _vel_3d);
    femData->outController.setOutput(var.name, var);
};

int FunctionPressureToElementVelocity::solveTimeStep(const NumLib::TimeStep &/*time*/)
{
    INFO("Calculating Darcy velocity within elements from fluid pressure...");
    getDarcyVelocity(*(NumLib::FemNodalFunctionScalar*)getInput(Pressure), *_vel);
    setOutput(Velocity, _vel_3d);
    return 0;
}

