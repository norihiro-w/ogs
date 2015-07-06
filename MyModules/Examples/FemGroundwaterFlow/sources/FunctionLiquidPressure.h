/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef LIQUID_PRESSURE_H
#define LIQUID_PRESSURE_H

#include <boost/property_tree/ptree.hpp>

#include "NumLib/Function/DiscreteDataConvergenceCheck.h"
#include "NumLib/Function/TXWrapped3DVectorFunction.h"
#include "NumLib/Function/FemIntegrationPointFunction.h"

#include "SolutionLib/Fem/FemIVBVProblem.h"
#include "SolutionLib/Fem/SingleStepFEM.h"

#include "THMCLib/ProcessLib/AbstractTransientProcess.h"


/**
 * \brief Liquid pressure calculator based on Groundwater flow equation using FEM
 *
 * Output result
 * - Liquid pressure (nodal)
 *
 * \tparam T_DISCRETE_SYSTEM    Discrete system type
 * \tparam T_LINEAR_SOLVER      Linear solver type
 */
class FunctionLiquidPressure
: public THMCLib::AbstractTransientProcess
{
public:
    // define input and output parameter id here
    enum Out { Pressure=0, Velocity = 1};

    ///
    FunctionLiquidPressure()
    : AbstractTransientProcess("LIQUID_FLOW", 0, 2), _problem(nullptr), _solution(nullptr), _feObjects(nullptr),
      _vel(nullptr), _vel_3d(nullptr)
    {
        // set default parameter name
        THMCLib::AbstractTransientProcess::setOutputParameterName(Pressure, "Pressure");
        THMCLib::AbstractTransientProcess::setOutputParameterName(Velocity, "Velocity");
    }

    ///
    virtual ~FunctionLiquidPressure()
    {
        delete _problem;
        delete _solution;
        delete _feObjects;
        delete _vel;
        delete _vel_3d;
    }

    /// initialize this process
    virtual bool initialize(const boost::property_tree::ptree&);

    /// finalize but nothing to do here
    virtual void finalize() {}

    /// return a convergence check class for this function
    virtual NumLib::IConvergenceCheck* getConvergenceChecker()
    {
        return &_checker;
    }

protected:
    virtual void postSolutionAlgorithm(const NumLib::TimeStep &/*time*/);

    virtual void updateOutputParameter(const NumLib::TimeStep &time);

    virtual SolutionLib::SingleStepFEM* getSolution() {return _solution;};

    virtual void output(const NumLib::TimeStep &time);


private:
    DISALLOW_COPY_AND_ASSIGN(FunctionLiquidPressure);

private:
    SolutionLib::FemIVBVProblem* _problem;
    SolutionLib::SingleStepFEM* _solution;
    NumLib::LagrangeFeObjectContainer* _feObjects;
    NumLib::DiscreteDataConvergenceCheck _checker;
    NumLib::FEMIntegrationPointFunctionVector* _vel;
    NumLib::TXWrapped3DVectorFunction<NumLib::FEMIntegrationPointFunctionVector>* _vel_3d;
};


#endif
