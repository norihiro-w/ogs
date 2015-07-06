/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 */

#pragma once

#include <boost/property_tree/ptree.hpp>

#include "NumLib/Function/DiscreteDataConvergenceCheck.h"
#include "NumLib/Function/TXWrapped3DVectorFunction.h"
#include "NumLib/Function/FemIntegrationPointFunction.h"

#include "SolutionLib/Fem/FemIVBVProblem.h"
#include "SolutionLib/Fem/SingleStepFEM.h"

#include "THMCLib/ProcessLib/AbstractTransientProcess.h"

/**
 *
 */
class FunctionTemperature
: public THMCLib::AbstractTransientProcess
{
public:
    enum In { Velocity = 0 };
    enum Out { Temperature = 0 };

    FunctionTemperature()
        : AbstractTransientProcess("HEAT_TRANSPORT", 1, 1),
          _problem(nullptr), _solution(nullptr), _feObjects(nullptr)
    {
        // set default parameter name
        THMCLib::AbstractTransientProcess::setInputParameterName(Velocity, "Velocity");
        THMCLib::AbstractTransientProcess::setOutputParameterName(Temperature, "Temperature");
    }

    virtual ~FunctionTemperature()
    {
        delete _problem;
        delete _solution;
        delete _feObjects;
    }

    /// initialize this process
    virtual bool initialize(const boost::property_tree::ptree &option);

    /// finalize but nothing to do here
    virtual void finalize() {};

    virtual NumLib::IConvergenceCheck* getConvergenceChecker()
    {
        return &_checker;
    }

protected:
    virtual void initializeTimeStep(const NumLib::TimeStep &time);

    virtual void updateOutputParameter(const NumLib::TimeStep &time);

    virtual SolutionLib::SingleStepFEM* getSolution() {return _solution;};

    virtual void output(const NumLib::TimeStep &time);

private:
    DISALLOW_COPY_AND_ASSIGN(FunctionTemperature);

private:
    SolutionLib::FemIVBVProblem* _problem;
    SolutionLib::SingleStepFEM* _solution;
    NumLib::LagrangeFeObjectContainer* _feObjects;
    NumLib::DiscreteDataConvergenceCheck _checker;
};

