################################################################################
# Configure subdirectory THMCLib
################################################################################

# Includes
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
  "${PROJECT_SOURCE_DIR}/THMCLib")
INCLUDE(SetupModules)

#-------------------------------------------------------------------------------
# Set up user modulues
#-------------------------------------------------------------------------------
SET(USER_MODULES_PATH ${CMAKE_SOURCE_DIR}/MyModules)
SETUP_USER_MODULES(${USER_MODULES_PATH} Module_SourceFiles TestModule_Files ALL_MODULE_TEST_NAME ALL_MODULE_TEST_DIR)
LIST(APPEND THMC_Files ${Module_SourceFiles}) 

#-------------------------------------------------------------------------------
# Set up other directories
#-------------------------------------------------------------------------------
# ./
FILE(GLOB TempProcessLib_Files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.h *.cpp *.txt *.cmake)
SOURCE_GROUP( ogs6THMC FILES ${TempProcessLib_Files} )
LIST(APPEND THMC_Files ${TempProcessLib_Files}) 

# ./modules
FILE(GLOB TempProcessLib_Files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} 
            ${USER_MODULES_PATH}/*.*)
SOURCE_GROUP( MyModules FILES ${TempProcessLib_Files} )
LIST(APPEND THMC_Files ${TempProcessLib_Files}) 

# ./MaterialLib
INCLUDE( MaterialLib/CMakeLists.txt )

# ./PhysicsLib
#INCLUDE( PhysicsLib/CMakeLists.txt )

# ./ChemLib
#INCLUDE( ChemLib/CMakeLists.txt )

# ./ProcessLib
INCLUDE( ProcessLib/CMakeLists.txt )

#
ADD_LIBRARY( THMCLib STATIC 
    ${MaterialLib_Files}
    ${ProcessLib_Files}
#    ${ChemLib_Files}
    ${PhysicsLib_Files}
    ${THMC_Files}
)

target_link_libraries( THMCLib
	BaseLib
	MathLib
	GeoLib
	MeshLib
	AssemblerLib
	NumLib
	SolutionLib
	FileIO
)

SET_TARGET_PROPERTIES(THMCLib PROPERTIES LINKER_LANGUAGE CXX)

INCLUDE_DIRECTORIES(
        .
        ${CMAKE_SOURCE_DIR} 
        ${USER_MODULES_PATH}
)
# This is for Configure.h which is generated later
INCLUDE_DIRECTORIES( ${OGS-Kernel_BINARY_DIR}/BaseLib )
# ProcessList.h and ProcessReg.h
INCLUDE_DIRECTORIES( ${PROJECT_BINARY_DIR} )

IF(USE_EXTERNAL_LIS)
    INCLUDE_DIRECTORIES( ${LIS_INCLUDE_DIR} )
ELSE()
	IF (USE_LIS)
		INCLUDE_DIRECTORIES( ${OGS-Kernel_SOURCE_DIR}/external/lis/include )
	ENDIF()
ENDIF ()

IF(USE_PETSC)
        TARGET_LINK_LIBRARIES( THMCLib  ${PETSC_LIBRARIES} )
    IF (MSVC)
        SET_TARGET_PROPERTIES(THMCLib PROPERTIES
                LINK_FLAGS "
            /NODEFAULTLIB:\"libcmt.lib\"
                ")
    ENDIF(MSVC)
ENDIF(USE_PETSC)

#-------------------------------------------------------------------------------
# Add tests
#-------------------------------------------------------------------------------
ADD_SUBDIRECTORY(Tests)

