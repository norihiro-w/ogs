/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 */

#include <gtest/gtest.h>
#include <logog/include/logog.hpp>
#include <logog/include/formatter.hpp>
#ifdef USE_LIS
#include <lis.h>
#endif
#ifdef USE_PETSC
#include <petscksp.h>
#endif

#include "TestEnvironment.h"

std::vector<std::string> g_args;

namespace
{
/**
 * new formatter for logog
 */
class FormatterCustom : public logog::FormatterGCC
{
    virtual TOPIC_FLAGS GetTopicFlags( const logog::Topic &topic )
    {
        return ( Formatter::GetTopicFlags( topic ) &
                 ~( TOPIC_LEVEL_FLAG | TOPIC_FILE_NAME_FLAG | TOPIC_LINE_NUMBER_FLAG ));
    }
};
}

int main(int argc, char *argv[])
{
    for (int i=0; i<argc; i++)
        g_args.push_back(argv[i]);

#ifdef USE_LIS
    lis_initialize(&argc, &argv);
#endif
#ifdef USE_PETSC
    char petsc_help[] = "Using PETSc package\n";
    PetscInitialize(&argc, &argv,(char *)0,petsc_help);
#endif

    int ret = 0;
    LOGOG_INITIALIZE();
//    try {
//        logog::Cout out;
//        FormatterCustom custom_format;
//        out.SetFormatter(custom_format);

        ::testing::InitGoogleTest(&argc, argv);
        ::testing::AddGlobalTestEnvironment(new TestEnvironment);
        ret = RUN_ALL_TESTS();
//    } catch (char* e) {
//        std::cerr << e << std::endl;
//    } catch (std::exception& e) {
//        std::cerr << e.what() << std::endl;
//    } catch (...) {
//        std::cerr << "Unknown exception occurred!" << std::endl;
//    }
    LOGOG_SHUTDOWN();

#ifdef USE_LIS
    lis_finalize();
#endif
#ifdef USE_PETSC
    PetscFinalize();
#endif

    return ret;
}

