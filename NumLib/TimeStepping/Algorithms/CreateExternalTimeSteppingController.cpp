/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "CreateExternalTimeSteppingController.h"
#include <string>

#include "BaseLib/ConfigTree.h"
#include "BaseLib/Error.h"
#include "BaseLib/FileTools.h"

#include "ExternalTimeSteppingController.h"
#include "TimeStepAlgorithm.h"

namespace NumLib
{
class TimeStepAlgorithm;
std::unique_ptr<TimeStepAlgorithm> createExternalTimeSteppingController(
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{prj__time_loop__processes__process__time_stepping__type}
    config.checkConfigParameter("type", "ExternalController");

    //! \ogs_file_param{prj__time_loop__processes__process__time_stepping__ExternalTimeSteppingController__t_initial}
    auto const t_initial = config.getConfigParameter<double>("t_initial");
    //! \ogs_file_param{prj__time_loop__processes__process__time_stepping__ExternalTimeSteppingController__t_end}
    auto const t_end = config.getConfigParameter<double>("t_end");
    //! \ogs_file_param{prj__time_loop__processes__process__time_stepping__ExternalTimeSteppingController__file}
    auto const filename = config.getConfigParameter<std::string>("file");
    auto const filepath = BaseLib::joinPaths(BaseLib::getProjectDirectory(), filename);

    auto opt_sleep = config.getConfigParameterOptional<unsigned>("sleep");
    unsigned const sleep = opt_sleep ? opt_sleep.get() : 100;
    if (sleep == 0)
        OGS_FATAL("Parameter <sleep> should be non-zero.");

    return std::make_unique<ExternalTimeSteppingController>(t_initial, t_end, filepath, sleep);
}
}  // end of namespace NumLib
