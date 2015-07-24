/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 */

#include "TimeSteppingController.h"

#include <iostream>

#include <logog/include/logog.hpp>

#include "BaseLib/MPITools.h"
#include "TimeStep.h"
#include "ITransientSystem.h"


namespace NumLib
{

size_t TimeSteppingController::solve(double time_end) 
{
    BaseLib::MPIEnvironment mpi;
    TimeStep time_current(_time_begin);

    while (time_current.current()<time_end) {
        double time_next = _root_subsystems->suggestNext(time_current);
        time_next = std::min(time_next, time_end);
        double suggested_dt = time_next - time_current.current();
        if (suggested_dt <= 0.0) {
            //error
            ERR("error - the suggested next time step is invalid.");
            break;
        }
        TimeStep t_n1(time_current);
        t_n1 += suggested_dt;

        mpi.barrier();
        if (mpi.root()) {
            INFO("\n");
            INFO("#############################################################");
            INFO("Time step %d: t=%f s, dt=%f s ", t_n1.steps(), time_next, t_n1.dt());
            INFO("#############################################################");
        }

        bool isAccepted = (_root_subsystems->solveTimeStep(t_n1)==0);
        if (isAccepted) {
            _root_subsystems->accept(t_n1);
            time_current += suggested_dt;
            doSomethingAfterTimeStepAccepted(time_current);
        }
    }

    return time_current.steps();
}

}
