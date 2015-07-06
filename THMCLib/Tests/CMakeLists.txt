################################################################################
# CMake script for testing modules
################################################################################

#-------------------------------------------------------------------------------
# Set up input and output path for testing
SET(USER_MODULES_TEST_OUTPUT_PATH "${PROJECT_BINARY_DIR}/TestMyModules")
# create directories for output
FOREACH (test_name ${ALL_MODULE_TEST_NAME})
    FILE(MAKE_DIRECTORY "${USER_MODULES_TEST_OUTPUT_PATH}/${test_name}")
ENDFOREACH()
# inform the root path to cpp
#ADD_DEFINITIONS("-DOGS_MY_MODULE_ROOT_DIR=${USER_MODULES_PATH}")
#ADD_DEFINITIONS("-DOGS_MY_MODULE_TEST_OUTPUT_ROOT_DIR=${USER_MODULES_TEST_OUTPUT_PATH}")
# create TestConfiguration.h
MESSAGE(STATUS "Creating TestConfiguration.h in ${PROJECT_BINARY_DIR}/MyModules")
SET(PROCESS_TEST_CONF_FILE ${PROJECT_BINARY_DIR}/MyModules/TestConfiguration.h)
FILE(WRITE ${PROCESS_TEST_CONF_FILE} "//Generated macros\n")
LIST(LENGTH ALL_MODULE_TEST_DIR len1)
MATH(EXPR len2 "${len1} - 1")

FOREACH(val RANGE ${len2})
    LIST(GET ALL_MODULE_TEST_NAME ${val} test_name)
    LIST(GET ALL_MODULE_TEST_DIR ${val} test_dir)
    # replace invalid characters for C macro and make it upper case
    STRING(REGEX REPLACE "/" "_" test_name_def ${test_name})
    STRING(REGEX REPLACE "-" "_" test_name_def ${test_name_def})
    STRING(TOUPPER ${test_name_def} test_name_def)
    # 
    FILE(APPEND ${PROCESS_TEST_CONF_FILE} "#define TEST_INPUT_PATH_${test_name_def} \"${test_dir}\"\n")
    FILE(APPEND ${PROCESS_TEST_CONF_FILE} "#define TEST_OUTPUT_PATH_${test_name_def} \"${USER_MODULES_TEST_OUTPUT_PATH}/${test_name}\"\n")
ENDFOREACH()

#-------------------------------------------------------------------------------
## pthread ##
SET ( CMAKE_THREAD_PREFER_PTHREAD On )
FIND_PACKAGE( Threads )
IF ( CMAKE_USE_PTHREADS_INIT )
        SET (HAVE_PTHREADS TRUE)
        MESSAGE (STATUS "pthread library found." )
ENDIF (CMAKE_USE_PTHREADS_INIT )

# source files
FILE(GLOB Tests_Files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.h *.cpp)
SOURCE_GROUP( testModules FILES ${Tests_Files} )

# Includes
INCLUDE_DIRECTORIES(
	.
	${PROJECT_BINARY_DIR}
)

# Create the executable
ADD_EXECUTABLE( test-ogs-THMC-Modules
        ${Tests_Files}
        ${TestModule_Files}
)

# Link
TARGET_LINK_LIBRARIES ( test-ogs-THMC-Modules
	BaseLib
	MathLib
	GeoLib
	MeshLib
	AssemblerLib
	NumLib
	SolutionLib
	THMCLib
    ${CMAKE_THREAD_LIBS_INIT}
    GTest
)

IF (WIN32)
        TARGET_LINK_LIBRARIES( test-ogs-THMC-Modules Winmm.lib )
ENDIF (WIN32)

IF(USE_OPENMP)
    FIND_PACKAGE( OpenMP REQUIRED )
    SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}" )
    SET( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}" )
    IF(NOT MSVC)
        SET( CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lgomp" )
    ENDIF(NOT MSVC)
ENDIF(USE_OPENMP)

# Properties
SET_TARGET_PROPERTIES(test-ogs-THMC-Modules PROPERTIES LINKER_LANGUAGE CXX)
