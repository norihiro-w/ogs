/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "ExternalTimeStepping.h"

#include <algorithm>
#include <cassert>
#include <limits>
#include <numeric>

#ifdef _WINDOWS
#include <windows.h>
#else
#include <unistd.h>
#define Sleep(x) usleep((x)*1000)
#endif

namespace NumLib
{
ExternalTimeStepping::ExternalTimeStepping(
    double t0, double tn, std::string const& timestep_file_path,
    unsigned sleep_duration_ms)
    : TimeStepAlgorithm(t0, tn),
      _timestep_file_path(timestep_file_path),
      _sleep_duration_ms(sleep_duration_ms)
{
}

bool ExternalTimeStepping::next(const double /*solution_error*/)
{
    // check if last time step
    if (std::abs(_ts_current.current() - _t_end) <
        std::numeric_limits<double>::epsilon())
        return false;

    // confirm current time and move to the next if accepted
    if (accepted())
    {
        _ts_prev = _ts_current;
        _dt_vector.push_back(_ts_current.dt());
    }

    // prepare the next time step info
    std::size_t time_step_index;
    double time_step_size;
    DBUG("Reading %s", _timestep_file_path.c_str());
    while (true)
    {
        std::ifstream ifs(_timestep_file_path);
        if (!ifs.good())
            OGS_FATAL("Failed to open %s", _timestep_file_path.c_str());

        ifs >> time_step_index >> time_step_size;

        if (time_step_index-1 == _ts_prev.steps())
            break;

        Sleep(_sleep_duration_ms);  // ms
    }

    if (std::abs(time_step_size) < std::numeric_limits<double>::epsilon())
        OGS_FATAL("Externally specified time step size is too small: dt=%g",
                  time_step_size);

    _ts_current = _ts_prev;
    _ts_current += time_step_size;

    return true;
}

void ExternalTimeStepping::finalizeCurrentTimeStep()
{
    std::ofstream ofs(_timestep_file_path);
    if (!ofs.good())
        OGS_FATAL("Failed to open %s", _timestep_file_path.c_str());

    ofs << _ts_current.steps() << "\n"
        << _ts_current.dt() << "\n"
        << 1 << std::endl;
}

}  // namespace NumLib
