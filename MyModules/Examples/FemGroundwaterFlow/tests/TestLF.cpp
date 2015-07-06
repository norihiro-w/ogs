/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 */

#include <string>

#include <gtest/gtest.h>

#include "BaseLib/FileTools.h"
#include "THMCLib/THMCSimulator.h"
#include "MyModules/TestConfiguration.h"
#include "ComparisonTools.h"
#include "ModuleTestTools.h"

#if defined(TEST_INPUT_PATH_FEMGROUNDWATERFLOW_LIQUID_FLOW_2D) \
    && defined(TEST_OUTPUT_PATH_FEMGROUNDWATERFLOW_LIQUID_FLOW_2D)

TEST(FemLiquidFlow, 2D)
{
    std::string str_in_dir = TEST_INPUT_PATH_FEMGROUNDWATERFLOW_LIQUID_FLOW_2D;
    std::string str_in_project = str_in_dir + "/q_quad";
    std::string str_out_dir = TEST_OUTPUT_PATH_FEMGROUNDWATERFLOW_LIQUID_FLOW_2D;
    std::vector<char*> argv(createCmdLineArguments(str_in_project, str_out_dir));
    {
        ogs6::THMCSimulator sim(argv.size(), &argv[0]);
        int ret = sim.execute();
        ASSERT_EQ(0, ret);
    }

    // compare results
    const std::string str_result_file = "q_quad_10.vtu";
    std::string str_ref_file = str_in_dir + "/reference_result/" + str_result_file;
    std::string str_new_file = str_out_dir + "/" + str_result_file;
    std::cout << "-> start comparing " << str_result_file << std::endl;
    bool ret_compare = compareAciiVTUFiles(str_ref_file, str_new_file, 1e-4);
    ASSERT_TRUE(ret_compare);
}

#endif
