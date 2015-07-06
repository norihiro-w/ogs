/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file SimulationInfo.cpp
 *
 * Created on 2012-07-05 by Norihiro Watanabe
 */


#pragma once

#include <string>
#include "BaseLib/CodingTools.h"

namespace ogs6
{

/**
 * \brief Simulation information
 */
class SimulationInfo
{
public:
    SimulationInfo();
    SimulationInfo(const std::string &project_path, const std::string &output_dir_path);

    static void outputHeader();

    void setProjectPath(const std::string &path);
    std::string getProjectPath() const {return _project_path;};
    std::string getProjectDirPath() const {return _project_dir;};
    std::string getProjectName() const {return _project_name;};
    std::string getOutputDirPath() const {return _output_dir;};

private:
    DISALLOW_COPY_AND_ASSIGN(SimulationInfo);

private:
    std::string _project_path;
    std::string _project_dir;
    std::string _project_name;
    std::string _output_dir;
};

}

