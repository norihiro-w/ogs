/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 */

#pragma once

#include <boost/property_tree/ptree.hpp>

#include "NumLib/Function/FemIntegrationPointFunction.h"
#include "NumLib/Function/DiscreteDataConvergenceCheck.h"
#include "NumLib/Function/TXWrapped3DVectorFunction.h"
#include "NumLib/TimeStepping/TimeStep.h"
#include "ProcessLib/AbstractTimeIndependentProcess.h"

class FunctionPressureToElementVelocity
    : public THMCLib::AbstractTimeIndependentProcess
{
public:
    enum In { Pressure=0 };
    enum Out { Velocity=0 };

    FunctionPressureToElementVelocity()
    : THMCLib::AbstractTimeIndependentProcess("PRESSURE_TO_ELEMENT_VELOCITY", 1, 1),
      _msh(nullptr), _vel(nullptr), _vel_3d(nullptr)
    {
        // set default parameter name
        THMCLib::AbstractTimeIndependentProcess::setInputParameterName(Pressure, "Pressure");
        THMCLib::AbstractTimeIndependentProcess::setOutputParameterName(Velocity, "Velocity");
    };

    virtual ~FunctionPressureToElementVelocity()
    {
        delete _vel;
        delete _vel_3d;
    };

    virtual bool initialize(const boost::property_tree::ptree &op);

    virtual void finalize() {};

    int solveTimeStep(const NumLib::TimeStep &/*time*/);

    virtual void accept(const NumLib::TimeStep &/*time*/);

    ///
    virtual NumLib::IConvergenceCheck* getConvergenceChecker() { return &_checker; };

private:
    MeshLib::Mesh* _msh;
    NumLib::FEMIntegrationPointFunctionVector* _vel;
    NumLib::DiscreteDataConvergenceCheck _checker;
    NumLib::TXWrapped3DVectorFunction<NumLib::FEMIntegrationPointFunctionVector>* _vel_3d;

    DISALLOW_COPY_AND_ASSIGN(FunctionPressureToElementVelocity);
};


