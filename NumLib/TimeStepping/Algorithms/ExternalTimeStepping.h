/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <memory>
#include <vector>

#include "TimeStepAlgorithm.h"

namespace NumLib
{

class ExternalTimeStepping final : public TimeStepAlgorithm
{
public:
    ExternalTimeStepping(
        double t0, double tn, std::string const& timestep_file_path,
        unsigned sleep_duration_ms);

    /// move to the next time step
    bool next(const double solution_error) override;

    /// return if current time step is accepted
    bool accepted() const override { return true; }

    void finalizeCurrentTimeStep() override;

private:
    std::string const _timestep_file_path;
    unsigned const _sleep_duration_ms; // in micro sec.
};

}  // NumLib
