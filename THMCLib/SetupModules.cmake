################################################################################
# CMake script: User modules parser
#
# This script assumes the following directory strucutre
#
# root/
# - MyModules/
#   - InactiveModules.cmake
#   - ModuleXXX/
#     - sources/
#     - tests/
#     - ConfigureModule.cmake
################################################################################

#-------------------------------------------------------------------------------
# Scans the current directory and returns a list of subdirectories.
# Third parameter is 1 if you want relative paths returned.
# Usage: list_subdirectories(the_list_is_returned_here C:/cwd 1)
#-------------------------------------------------------------------------------
MACRO(LIST_SUBDIRECTORIES retval curdir return_relative)
#    MESSAGE (STATUS "curdir: ${curdir}")
#    file(GLOB sub-dir ${curdir}/*)
    FILE(GLOB sub-dir RELATIVE ${curdir} ${curdir}/*)
#    list(LENGTH sub-dir n_dir)
#    MESSAGE (STATUS "sub-dir (${n_dir}): ${sub-dir}")
    SET(list_of_dirs "")
    FOREACH(dir ${sub-dir})
        IF(IS_DIRECTORY ${curdir}/${dir})
#            MESSAGE (STATUS "${dir} is a directory")
            IF (${return_relative})
                SET(list_of_dirs ${list_of_dirs} ${dir})
            ELSE()
                SET(list_of_dirs ${list_of_dirs} ${curdir}/${dir})
            ENDIF()
        ENDIF()
    ENDFOREACH()
    SET(${retval} ${list_of_dirs})
ENDMACRO()

#-------------------------------------------------------------------------------
# Add prefix to each entriy in a list
#-------------------------------------------------------------------------------
MACRO(LIST_ADD_PREFIX retval curlist new_prefix)
    SET(new_list)
    FOREACH(dir ${${curlist}})
        SET(new_list ${new_list} ${new_prefix}${dir})
    ENDFOREACH()
    SET(${retval} ${new_list})
ENDMACRO()

#-------------------------------------------------------------------------------
# Add macro descriptions to register processes
# - OGS_ADD_PROCESS_SYS_SOLVER
# - OGS_ADD_PROCESS_SYS
# - OGS_ADD_PROCESS
#-------------------------------------------------------------------------------
MACRO(OGS_ADD_PROCESS_SYS_SOLVER 
    PROCESS_LIST        # Process list (output)
    PROCESS_NAME        # Process name
    PROCESS_CLASS       # Complete class name with namespace 
    DISCRETE_SYS_TYPE   # Discrete system type
    LINEAR_SOLVER_TYPE  # Linear solver type
    )
    SET(TMP_PROCESS_REGIST 
        "OGS_ADD_PROCESS_SYS_SOLVER(${PROCESS_NAME}, ${PROCESS_CLASS}, ${DISCRETE_SYS_TYPE}, ${LINEAR_SOLVER_TYPE});"
    )
#    MESSAGE (STATUS "TMP_PROCESS_REGIST is ${TMP_PROCESS_REGIST}")
    SET(${PROCESS_LIST} ${${PROCESS_LIST}} ${TMP_PROCESS_REGIST})
#    MESSAGE (STATUS "OGS_PROCESS_LIST is ${${PROCESS_LIST}}")
ENDMACRO()

MACRO(OGS_ADD_PROCESS_SYS 
    PROCESS_LIST        # Process list (output)
    PROCESS_NAME        # Process name
    PROCESS_CLASS       # Complete class name with namespace 
#    DISCRETE_SYS_TYPE   # Discrete system type
    )
    SET(TMP_PROCESS_REGIST 
        "OGS_ADD_PROCESS_SYS(${PROCESS_NAME}, ${PROCESS_CLASS});"
#        "OGS_ADD_PROCESS_SYS(${PROCESS_NAME}, ${PROCESS_CLASS}, ${DISCRETE_SYS_TYPE});"
    )
    SET(${PROCESS_LIST} ${${PROCESS_LIST}} ${TMP_PROCESS_REGIST})
ENDMACRO()

MACRO(OGS_ADD_PROCESS 
    PROCESS_LIST        # Process list (output)
    PROCESS_NAME        # Process name
    PROCESS_CLASS       # Complete class name with namespace 
    )
    SET(TMP_PROCESS_REGIST 
        "OGS_ADD_PROCESS(${PROCESS_NAME}, ${PROCESS_CLASS});"
    )
    SET(${PROCESS_LIST} ${${PROCESS_LIST}} ${TMP_PROCESS_REGIST})
ENDMACRO()

#-------------------------------------------------------------------------------
# Parse one module
#
# input:
# - MODULE_DIR
# - MODULE_NAME
#
# output:
# - VALID_MODULE
# - MODULE_PROCESS_HEADERS
# - MODULE_PROCESS_LIST
# - MODULE_SOURCE_FILES
# - MODULE_TEST_FILES
# - MODULE_TEST_NAME
# - MODULE_TEST_INPUT_DIR
#-------------------------------------------------------------------------------
MACRO(PARSE_MODULE MODULE_DIR MODULE_NAME VALID_MODULE MODULE_PROCESS_HEADERS MODULE_PROCESS_LIST MODULE_SOURCE_FILES MODULE_TEST_FILES MODULE_TEST_NAME MODULE_TEST_INPUT_DIR)

    #MESSAGE(STATUS "Module: ${MODULE_DIR}")

    # need this. otherwise this variable is not paased to ConfigureModule.cmake
    SET(MODULE_DIR ${MODULE_DIR})
    SET(MODULE_NAME ${MODULE_NAME})
    
    # include to check MODULE_SUBDIRECTORIES
    SET(MODULE_SUBDIRECTORIES) #reset
    INCLUDE(${MODULE_DIR}/ConfigureModule.cmake)
    
    IF (MODULE_SUBDIRECTORIES) # if this directory has modules in subdirectories
        APPEND_MODULE_DIRECTORIES(${MODULE_DIR} MODULE_PROCESS_HEADERS MODULE_PROCESS_LIST MODULE_SOURCE_FILES MODULE_TEST_FILES MODULE_TEST_NAME MODULE_TEST_INPUT_DIR)
        SET(VALID_MODULE FALSE)
    ELSE()
        # reset var
        SET(MODULE_PROCESS_HEADERS)
        SET(MODULE_PROCESS_LIST)
        SET(MODULE_SOURCE_FILES)
        SET(MODULE_TEST_FILES)
        SET(MODULE_TEST_INPUT_DIR)
        SET(MODULE_TEST_NAME)
        SET(MODULE_SUBDIRECTORIES)
        
        # check this module
        SET(VALID_MODULE TRUE)
        IF(NOT EXISTS ${MODULE_DIR}/ConfigureModule.cmake)
            MESSAGE(STATUS "*** Warning: ConfigureModule.cmake is not found in a module ${dir}. This module is not loaded.")
            SET(VALID_MODULE FALSE)
        ENDIF()
        IF(NOT EXISTS ${MODULE_DIR}/sources)
            MESSAGE(STATUS "*** Warning: a subdirectory 'sources' is not found in a module ${dir}. This module is not loaded.")
            SET(VALID_MODULE FALSE)
        ENDIF()
        IF(NOT EXISTS ${MODULE_DIR}/tests)
            MESSAGE(STATUS "*** Warning: a subdirectory 'tests' is not found in a module ${dir}. ")
        ENDIF()
        
        # load this module     
        IF(VALID_MODULE)
            # include
            INCLUDE(${MODULE_DIR}/ConfigureModule.cmake)
            IF(NOT MODULE_SOURCE_FILES)
				MESSAGE(STATUS "${MODULE_DIR}")
                # get module source files under sources/            
				SET(MODULE_SRC_DIR "${MODULE_DIR}/sources")
                FILE(GLOB MODULE_SOURCE_FILES "${MODULE_SRC_DIR}/*")
				#MESSAGE(STATUS "- ${MODULE_SOURCE_FILES}")
				LIST(LENGTH ${MODULE_SOURCE_FILES} len)
				math(EXPR _stop "${len}-1")
				MESSAGE(STATUS "- ${len} source files found")
				#MESSAGE(STATUS "- ${inFile}")
				FOREACH(idx RANGE 0 ${_stop})
					LIST(GET ${MODULE_SOURCE_FILES} ${idx} inFile)
					#MESSAGE(STATUS "- ${idx}: ${inFile}")
					SOURCE_GROUP( ${SOURCE_GROUP_USER_MODULES_SOURCES_ROOT}\\${MODULE_NAME} FILES ${inFile} )
				ENDFOREACH()
                #SOURCE_GROUP( ${SOURCE_GROUP_USER_MODULES_SOURCES_ROOT}\\${MODULE_NAME} FILES ${MODULE_SOURCE_FILES} )
            ENDIF()
            IF(NOT MODULE_TEST_FILES)
                # get module test source files under tests/            
                FILE(GLOB MODULE_TEST_FILES ${MODULE_DIR}/tests/*.*)
                #SOURCE_GROUP( ${SOURCE_GROUP_USER_MODULES_TEST_ROOT}\\${MODULE_NAME} FILES ${MODULE_TEST_FILES} )
				LIST(LENGTH ${MODULE_TEST_FILES} len)
				math(EXPR _stop "${len}-1")
				MESSAGE(STATUS "- ${len} test files found")
				FOREACH(idx RANGE 0 ${_stop})
					LIST(GET ${MODULE_TEST_FILES} ${idx} inFile)
					#MESSAGE(STATUS "- ${idx}: ${inFile}")
					SOURCE_GROUP( ${SOURCE_GROUP_USER_MODULES_TEST_ROOT}\\${MODULE_NAME} FILES ${inFile} )
				ENDFOREACH()
				#MESSAGE(STATUS "${SOURCE_GROUP_USER_MODULES_TEST_ROOT}\\${MODULE_NAME}")
            ENDIF()
            #
            #MESSAGE(STATUS "MODULE_PROCESS_HEADERS=${${MODULE_PROCESS_HEADERS}}")
            LIST_ADD_PREFIX(MODULE_PROCESS_HEADERS ${MODULE_PROCESS_HEADERS} "${MODULE_DIR}/sources/")
            LIST_ADD_PREFIX(MODULE_TEST_NAME ${MODULE_TEST_INPUT_DIR} "${MODULE_NAME}/")
            LIST_ADD_PREFIX(MODULE_TEST_INPUT_DIR ${MODULE_TEST_INPUT_DIR} "${MODULE_DIR}/tests/")
        ENDIF()
    ENDIF()


ENDMACRO()

#-------------------------------------------------------------------------------
# Parse a directory with modules 
#
# input:
# - ROOT_DIR
#
# output:
# - ALL_MODULE_HEADER_FILES 
# - ALL_OGS_PROCESS_LIST 
# - THMC_Files 
# - TestModule_Files 
# - ALL_MODULE_TEST_NAME
# - ALL_MODULE_TEST_DIR
#-------------------------------------------------------------------------------
MACRO(APPEND_MODULE_DIRECTORIES ROOT_DIR ALL_MODULE_HEADER_FILES ALL_OGS_PROCESS_LIST THMC_Files TestModule_Files ALL_MODULE_TEST_NAME ALL_MODULE_TEST_DIR)
    MESSAGE (STATUS "Searching module directories in ${ROOT_DIR}")
    LIST_SUBDIRECTORIES(List_subdir ${ROOT_DIR} 1)
    LIST(LENGTH List_subdir n_dir)
    MESSAGE (STATUS "${n_dir} directories are found.")
    
    # Exclude modules given in InactiveModules.cmake 
    IF(EXISTS ${ROOT_DIR}/InactiveModules.cmake)
        INCLUDE(${ROOT_DIR}/InactiveModules.cmake)
        LIST(LENGTH INACTIVE_MODULES n_inactive_dir)
        IF (n_inactive_dir GREATER 0)
            MESSAGE (STATUS "${n_inactive_dir} modules are excluded. Inactive modules are defined in ${ROOT_DIR}/InactiveModules.cmake")
            LIST(REMOVE_ITEM List_subdir ${INACTIVE_MODULES})
        ENDIF()
    ENDIF()
    
    # Configure active modules
    FOREACH(dir ${List_subdir})
    #    MESSAGE (STATUS "sub dir found: ${dir}")
        PARSE_MODULE(${ROOT_DIR}/${dir} ${dir} VALID_MODULE MODULE_PROCESS_HEADERS MODULE_PROCESS_LIST MODULE_SOURCE_FILES MODULE_TEST_FILES MODULE_TEST_NAME MODULE_TEST_INPUT_DIR)
        IF (VALID_MODULE)
            # Add process 
            LIST(APPEND ALL_MODULE_HEADER_FILES ${MODULE_PROCESS_HEADERS})
            LIST(APPEND ALL_OGS_PROCESS_LIST ${MODULE_PROCESS_LIST})
            LIST(APPEND THMC_Files ${MODULE_SOURCE_FILES}) 
            # Add test
            LIST(APPEND TestModule_Files ${MODULE_TEST_FILES})
            LIST(APPEND ALL_MODULE_TEST_NAME ${MODULE_TEST_NAME})
            LIST(APPEND ALL_MODULE_TEST_DIR ${MODULE_TEST_INPUT_DIR})
        ENDIF()
#        LIST(LENGTH ALL_OGS_PROCESS_LIST n_pcs)
#        MESSAGE (STATUS "${n_pcs} pcs are found so far. ${dir}")
    ENDFOREACH()
ENDMACRO()

#-------------------------------------------------------------------------------
# Set up user modules
# 
# This macro parses all directories under MyModules to find user-defined modules  
# and return some variables. It also generates ProcessList.h and ProcessReg.h
# in ${PROJECT_BINARY_DIR}/MyModules.
# 
# input:
# - USER_MODULES_PATH
#
# output:
# - THMC_Files
# - TestModule_Files
# - ALL_MODULE_TEST_NAME
# - ALL_MODULE_TEST_DIR
#-------------------------------------------------------------------------------
MACRO(SETUP_USER_MODULES USER_MODULES_PATH THMC_Files TestModule_Files ALL_MODULE_TEST_NAME  ALL_MODULE_TEST_DIR)

# init
SET(SOURCE_GROUP_USER_MODULES_SOURCES_ROOT MyModules)
SET(SOURCE_GROUP_USER_MODULES_TEST_ROOT testModules)

# Look for directories (i.e. modules)
APPEND_MODULE_DIRECTORIES(${USER_MODULES_PATH} ALL_MODULE_HEADER_FILES ALL_OGS_PROCESS_LIST THMC_Files TestModule_Files ALL_MODULE_TEST_NAME ALL_MODULE_TEST_DIR)

# Exclude modules given in InactiveModules.cmake 
IF(EXISTS ${USER_MODULES_PATH}/InactiveModules.cmake)
    INCLUDE(${USER_MODULES_PATH}/InactiveModules.cmake)
    LIST(LENGTH INACTIVE_MODULES n_inactive_dir)
    IF (n_inactive_dir GREATER 0)
        MESSAGE (STATUS "${n_inactive_dir} modules are excluded. Inactive modules are defined in MyModules/InactiveModules.cmake")
        LIST(REMOVE_ITEM List_subdir ${INACTIVE_MODULES})
    ENDIF()
ENDIF()

# for debugging
LIST(LENGTH ALL_OGS_PROCESS_LIST n_pcs)
MESSAGE (STATUS "${n_pcs} processes will be registered.")
#MESSAGE (STATUS "${ALL_OGS_PROCESS_LIST}")

# Create a process header file
MESSAGE(STATUS "Creating ProcessList.h and ProcessReg.h in ${PROJECT_BINARY_DIR}/MyModules")
SET(PROCESS_LIST_FILE ${PROJECT_BINARY_DIR}/MyModules/ProcessList.h)
FILE(WRITE ${PROCESS_LIST_FILE} "//Generated headers include\n")
FOREACH(header ${ALL_MODULE_HEADER_FILES})
    FILE(APPEND ${PROCESS_LIST_FILE} "#include \"${header}\"\n")
ENDFOREACH()

# Create a process registration file
SET(PROCESS_REGISTER_FILE ${PROJECT_BINARY_DIR}/MyModules/ProcessReg.h)
#MESSAGE(STATUS "Creating ProcessReg.h in ${PROCESS_REGISTER_FILE}")
FILE(WRITE ${PROCESS_REGISTER_FILE} "//Generated process registration\n")
FOREACH(process_regist ${ALL_OGS_PROCESS_LIST})
    FILE(APPEND ${PROCESS_REGISTER_FILE} "${process_regist};\n")
ENDFOREACH()

ENDMACRO()
