#####################
### OGS-6 Project ###
#####################

# Specify minimum CMake version
cmake_minimum_required(VERSION 2.6)

# Project name
project( OGS-6 )

### CMake includes ###
INCLUDE(scripts/cmake/CheckTypeSizes.cmake)
INCLUDE(scripts/cmake/FindIncludeHeader.cmake)
INCLUDE(scripts/cmake/Functions.cmake)
INCLUDE(scripts/cmake/CMakeSetup.cmake)
INCLUDE(scripts/cmake/CompilerSetup.cmake)
INCLUDE(scripts/cmake/Find.cmake)
INCLUDE(scripts/cmake/ProjectSetup.cmake)
INCLUDE(scripts/cmake/DocumentationSetup.cmake)

###############
### Options ###
###############

# Configurations
OPTION(USE_OPENMP "Use OpenMP" OFF)
OPTION(USE_MPI "Use MPI" OFF)
#OPTION(USE_EIGEN "Use Eigen" OFF)
OPTION(USE_LIS "Use LIS" ON)
OPTION(USE_EXTERNAL_LIS "Use external LIS" OFF)
OPTION(USE_BLAS_LAPACK "Use BLAS&LAPACK" OFF)
OPTION(USE_PETSC "Use PETSc libs for parallel computing" OFF)
OPTION(USE_MKL_PARDISO "Use PARDISO in MKL library" OFF)
OPTION(USE_ORIGINAL_PARDISO "Use original PARDISO solver" OFF)
OPTION(PYTHON_WRAP "Provide Python interface" OFF)
OPTION(OGS_BUILD_INFO "Should build information be generated" ON)
OPTION(BUILD_OGS_TEST "Should build test" OFF)

IF(USE_MKL_PARDISO AND USE_ORIGINAL_PARDISO)
    MESSAGE (FATAL_ERROR "Aborting: Please choose USE_MKL_PARDISO or USE_ORIGINAL_PARDISO. Selecting both is not supported.")
ENDIF()

###########################################################################
### OGS version information. Adjust these if you release a new version. ###
###########################################################################
SET (OGS_VERSION_MAJOR 0)
SET (OGS_VERSION_MINOR 1)
SET (OGS_VERSION_PATCH 0)
SET (OGS_RELEASE_PERSONS "NW")
SET (OGS_VERSION "${OGS_VERSION_MAJOR}.${OGS_VERSION_MINOR}.${OGS_VERSION_PATCH}(${OGS_RELEASE_PERSONS})")
# should we use ISO dates: 2010-12-31 ?
SET (OGS_DATE "05.07.2012")


###########################################################################

IF(USE_OPENMP)
	FIND_PACKAGE( OpenMP REQUIRED )
	SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}" )
	SET( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}" )
	IF(NOT MSVC)
		SET( CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lgomp" )
	ENDIF(NOT MSVC)
ENDIF(USE_OPENMP)

IF (USE_MPI)
	FIND_PACKAGE(MPI)
	IF(MPI_FOUND)		
	  SET(CMAKE_C_COMPILER ${MPI_COMPILER})
	  SET(CMAKE_CXX_COMPILER ${MPI_COMPILER})
		ADD_DEFINITIONS(-DUSE_MPI)
	ELSE(MPI_FOUND)
	  MESSAGE (FATAL_ERROR "Aborting: MPI implementation is not found!")
	ENDIF(MPI_FOUND)			
ENDIF (USE_MPI)

INCLUDE_DIRECTORIES( external )
#IF(USE_EIGEN)
#	FIND_PACKAGE ( Eigen REQUIRED )
	ADD_DEFINITIONS(-DUSE_EIGEN)
#ENDIF(USE_EIGEN)

IF(USE_LIS)
	IF (USE_EXTERNAL_LIS)
		FIND_PACKAGE ( LIS REQUIRED )
	ENDIF()
	ADD_DEFINITIONS(-DUSE_LIS)
ENDIF(USE_LIS)

IF(USE_BLAS_LAPACK)
	ADD_DEFINITIONS(-DUSE_BLAS_LAPACK)
	IF(NOT MSVC)
		FIND_PACKAGE ( BLAS REQUIRED )
		FIND_PACKAGE ( LAPACK REQUIRED )
	ELSE(NOT MSVC)
		find_library(BLAS_LIBRARIES
			NAMES libf2cblas.lib
			PATHS ${CMAKE_SOURCE_DIR}/../ )	
		find_library(LAPACK_LIBRARIES
			NAMES libf2clapack.lib
			PATHS ${CMAKE_SOURCE_DIR}/../ )	
	ENDIF(NOT MSVC)

ENDIF(USE_BLAS_LAPACK)

IF(USE_PETSC)
	MESSAGE (STATUS  "Configuring for PETSc" )
	ADD_DEFINITIONS(-DUSE_PETSC)
#	IF(UNIX)
#		ADD_DEFINITIONS(-Wunused-parameter)
#	ENDIF(UNIX)
	FIND_PACKAGE(PETSc REQUIRED)
	IF (NOT USE_MPI)
		set (PETSC_INCLUDES ${PETSC_DIR}/include/mpiuni ${PETSC_INCLUDES})
	ENDIF (NOT USE_MPI)
	INCLUDE_DIRECTORIES( ${PETSC_INCLUDES} )
ENDIF(USE_PETSC)

IF(USE_MKL_PARDISO)
    SET(USE_MKL ON)
    SET(USE_PARDISO ON)
    ADD_DEFINITIONS(-DUSE_MKL_PARDISO)
ENDIF(USE_MKL_PARDISO)

IF(USE_ORIGINAL_PARDISO)
    SET(USE_PARDISO ON)
    ADD_DEFINITIONS(-DUSE_ORIGINAL_PARDISO)
ENDIF(USE_ORIGINAL_PARDISO)

IF (USE_PARDISO)
    ADD_DEFINITIONS(-DUSE_PARDISO)
ENDIF (USE_PARDISO)

IF(USE_MKL)
    MESSAGE (STATUS  "Configuring for MKL" )
    FIND_PACKAGE(MKL REQUIRED)
    INCLUDE_DIRECTORIES( ${MKL_INCLUDES} )
ENDIF(USE_MKL)

IF (PYTHON_WRAP)
	FIND_PACKAGE(SWIG REQUIRED)
	INCLUDE(${SWIG_USE_FILE})
	FIND_PACKAGE(PythonLibs)
	INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})
ENDIF()

# Profiling
IF((CMAKE_COMPILER_IS_GNUCXX OR CMAKE_COMPILER_IS_GNUCC) AND GPROF_PATH)
	OPTION(OGS_PROFILE "Enables compiling with flags set for profiling with gprof." OFF)
ENDIF() # GCC AND GPROF_PATH


# Print CMake variable values
IF (OGS_CMAKE_DEBUG)
	INCLUDE(ListAllCMakeVariableValues)
	list_all_cmake_variable_values()
ENDIF (OGS_CMAKE_DEBUG)

######################
### Subdirectories ###
######################

# Add subdirectories with the projects
ADD_SUBDIRECTORY( ogsCore )
ADD_SUBDIRECTORY( ogsNumerics )
ADD_SUBDIRECTORY( ogsFileIO )
ADD_SUBDIRECTORY( ogs6THMC )
ADD_SUBDIRECTORY( ogs6CUI )
ADD_SUBDIRECTORY( external )
IF (BUILD_OGS_TEST)
    IF (NOT GTEST_FOUND)
        find_package(GTest)
    ENDIF()
    IF (NOT GTEST_FOUND)
          MESSAGE (FATAL_ERROR "Aborting: Gtest lib is not found!")
    ENDIF()    
	INCLUDE_DIRECTORIES(${GTEST_INCLUDE_DIRS})
	ADD_SUBDIRECTORY(Tests)
ENDIF ()
IF(USE_MPI)
	ADD_SUBDIRECTORY(TestsMPI)
ENDIF(USE_MPI)
IF (PYTHON_WRAP)
	ADD_SUBDIRECTORY(ogs6PythonUI)
ENDIF()

# Create Configure.h
CONFIGURE_FILE (ogsCore/BaseLib/Configure.h.in ${PROJECT_BINARY_DIR}/ogsCore/BaseLib/Configure.h)
CONFIGURE_FILE (ogsCore/BaseLib/BuildInfo.h.in ${PROJECT_BINARY_DIR}/ogsCore/BaseLib/BuildInfo.h)
