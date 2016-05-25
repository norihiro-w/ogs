#####################
### OGS-6 Project ###
#####################

# Specify minimum CMake version
cmake_minimum_required(VERSION 3.1)

# Set CMake policies
cmake_policy(SET CMP0011 OLD)
cmake_policy(SET CMP0054 NEW)

# Project name
project( OGS-6 )



#TODO
#- set OGS_USE_XXX=ON if a package is found
#- use REQUIRED to FindPackage if OGS_USE_XXX=ON
#- cmake file for each module? OGS_USE_XXX, find, compiler options, ?

####################################
### OGS options                  ###
####################################

#------------------------------------
# Build targets
#------------------------------------
option(OGS_BUILD_CLI "Should the OGS simulator be built?" ON)
option(OGS_BUILD_GUI "Should the Data Explorer be built?" OFF)
option(OGS_BUILD_TESTS "Should the test executables be built?" ON)
option(OGS_BUILD_UTILS "Should the utilities programms be built?" OFF)

#------------------------------------
# Compiling options
#------------------------------------
# Compiler flags
set(OGS_CXX_FLAGS "" CACHE STRING "Additional C++ compiler flags.")
option(STL_NO_DEBUG "Disable STL debug in debug build" OFF)
set(OGS_CPU_ARCHITECTURE "native" CACHE STRING "Processor architecture, defaults to native.")

# Logging
option(OGS_DISABLE_LOGGING "Disables all logog messages." OFF)

# Paralleization
option(OGS_USE_MPI "Use MPI" OFF)
option(OGS_USE_OPENMP "Use OpenMP" OFF)

# Options controlling which FEM elements will be compiled
set(OGS_MAX_ELEMENT_DIM   3 CACHE STRING "Maximum dimension of FEM elements to be built.")
set(OGS_MAX_ELEMENT_ORDER 2 CACHE STRING "Maximum order of FEM elements to be built.")
option(OGS_ENABLE_ELEMENT_SIMPLEX "Build FEM elements for simplices (triangles, tetrahedra)." ON)
option(OGS_ENABLE_ELEMENT_CUBOID  "Build FEM elements for cuboids (quads, hexahedra)." ON)
option(OGS_ENABLE_ELEMENT_PRISM   "Build FEM elements for prisms." ON)
option(OGS_ENABLE_ELEMENT_PYRAMID "Build FEM elements for pyramids." ON)

#------------------------------------
# External libraries
#------------------------------------

option(OGS_NO_EXTERNAL_LIBS "Builds OGS without any external dependencies." OFF)

option(OGS_USE_CVODE "Use CVODE" OFF)
option(OGS_USE_LAPACK "Use LAPACL" OFF)
option(OGS_USE_METIS "Use METIS" OFF)
option(OGS_USE_SHAPELIB "Use ShapeLib" OFF)
option(OGS_USE_QT "Use Qt" OFF)

# Eigen
option(OGS_USE_EIGEN "Use EIGEN for local matrix and vector" ON)
option(OGS_EIGEN_DYNAMIC_SHAPE_MATRICES "Use dynamically allocated shape matrices" ON)
option(EIGEN_NO_DEBUG "Disables Eigen's assertions" OFF)

# Linear solvers
option(OGS_USE_LIS "Use Lis" OFF)
option(OGS_USE_EIGENLIS "Use Lis solver together with Eigen sparse matrices" OFF)
option(OGS_USE_PETSC "Use PETSc routines" OFF)

# Third-party libraries
# Default checks for system first, then builds locally
set(OGS_LIB_VTK "Default" CACHE STRING "Which VTK library should be used?")
set(OGS_LIB_BOOST "Default" CACHE STRING "Which Boost library should be used?")
set(OGS_LIB_EIGEN "Default" CACHE STRING "Which Eigen library should be used?")
set_property(CACHE
    OGS_LIB_VTK
    OGS_LIB_BOOST
    OGS_LIB_EIGEN
    PROPERTY STRINGS "Default" "System" "Local")


#------------------------------------
# profiling, packaging, etc.
#------------------------------------

# Code coverage
option(OGS_COVERAGE "Enables code coverage measurements with gcov/lcov." OFF)

# Profiling
option(OGS_PROFILE "Enables compiling with flags set for profiling with gprof." OFF)

# Packaging
option(OGS_DOWNLOAD_ADDITIONAL_CONTENT "Should addional content such as manuals be downloaded and packaged?" OFF)


#######################################
### OGS options overwrite           ###
#######################################

if(OGS_BUILD_GUI)
    set(OGS_USE_QT ON CACHE BOOL "Use Qt" FORCE)
    set(OGS_USE_SHAPELIB ON CACHE BOOL "Use ShapeLib" FORCE)
endif()

if (OGS_USE_EIGENLIS)
    set(OGS_USE_EIGEN ON CACHE BOOL "Use Eigen" FORCE)
    set(OGS_USE_LIS ON CACHE BOOL "Use Lis" FORCE)
endif()


###########################
### Setup               ###
###########################
set(CMAKE_LIBRARY_SEARCH_PATH "" CACHE PATH
    "Additional library installation path, e.g. /opt/local or C:/libs")
list(APPEND CMAKE_PREFIX_PATH
    $ENV{HOMEBREW_ROOT}             # Homebrew package manager on Mac OS
    $ENV{CMAKE_LIBRARY_SEARCH_PATH} # Environment variable, Windows
    ${CMAKE_LIBRARY_SEARCH_PATH})   # CMake option, Windows

### CMake includes ###
include(scripts/cmake/ConanSetup.cmake)
include(scripts/cmake/CheckTypeSizes.cmake)
include(scripts/cmake/Functions.cmake)
include(scripts/cmake/CMakeSetup.cmake)
include(scripts/cmake/CompilerSetup.cmake)
include(scripts/cmake/Find.cmake)
include(scripts/cmake/SubmoduleSetup.cmake)
include(scripts/cmake/ProjectSetup.cmake)
include(scripts/cmake/DocumentationSetup.cmake)
include(scripts/cmake/test/Test.cmake)
if(OGS_COVERAGE AND NOT IS_SUBPROJECT)
    include(scripts/cmake/Coverage.cmake)
endif()

# Print CMake variable values
if (OGS_CMAKE_DEBUG)
    include(ListAllCMakeVariableValues)
    list_all_cmake_variable_values()
endif ()

include(scripts/cmake/packaging/Pack.cmake)


###################
### Definitions ###
###################
if(OGS_BUILD_GUI)
    add_definitions(-DOGS_BUILD_GUI)
endif()

if(OGS_BUILD_TESTS)
    set(Data_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/Tests/Data CACHE INTERNAL "")
    set(Data_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/Tests/Data CACHE INTERNAL "")
endif()

if(OGS_USE_CVODE)
    add_definitions(-DCVODE_FOUND)
endif()

if(OGS_USE_LIS)
    add_definitions(-DUSE_LIS)
    include_directories(SYSTEM ${LIS_INCLUDE_DIR})
endif()

if(OGS_USE_PETSC)
    add_definitions(-DUSE_PETSC)
    set(OGS_USE_MPI ON CACHE BOOL "Use MPI" FORCE)
endif()

if(OGS_USE_MPI)
    add_definitions(-DUSE_MPI)
    include_directories(SYSTEM ${MPI_CXX_INCLUDE_PATH})
endif()

if(OGS_USE_OPENMP)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    message(STATUS "OpenMP enabled.")
endif()

if(OGS_USE_EIGEN)
    add_definitions(-DEIGEN_INITIALIZE_MATRICES_BY_ZERO) # TODO check if needed
    add_definitions(-DOGS_USE_EIGEN)
    if (EIGEN_NO_DEBUG)
        add_definitions(-DEIGEN_NO_DEBUG)
    endif()
    if(OGS_EIGEN_DYNAMIC_SHAPE_MATRICES)
        add_definitions(-DOGS_EIGEN_DYNAMIC_SHAPE_MATRICES)
    endif()
endif()

if (OGS_USE_EIGENLIS)
    add_definitions(-DOGS_USE_EIGENLIS)
endif()

if(OGS_USE_SHAPELIB)
    include_directories(SYSTEM ${Shapelib_INCLUDE_DIRS})
endif()


######################
### Subdirectories ###
######################

# External projects
include(scripts/cmake/ExternalProjectBoost.cmake)
include(VtkModules)
if(OGS_INSITU)
    include(ExternalProjectCatalyst)
else()
    include(ExternalProjectVtk)
endif()

include_directories( SYSTEM ${Boost_INCLUDE_DIRS} )
if(OGS_USE_EIGEN)
    include(scripts/cmake/ExternalProjectEigen.cmake)
    include_directories(SYSTEM ${EIGEN3_INCLUDE_DIR})
endif()

# Add subdirectories with the projects
add_subdirectory( ThirdParty )
include_directories( SYSTEM ${CMAKE_CURRENT_SOURCE_DIR}/ThirdParty )
include_directories( SYSTEM ${CMAKE_CURRENT_SOURCE_DIR}/ThirdParty/gtest/include )
include_directories( SYSTEM ${CMAKE_CURRENT_SOURCE_DIR}/ThirdParty/autocheck/include )
include_directories( SYSTEM ${CMAKE_CURRENT_SOURCE_DIR}/ThirdParty/tclap/include )
include_directories(SYSTEM ${CMAKE_CURRENT_SOURCE_DIR}/ThirdParty/vtkGUISupportQt)

add_subdirectory( Applications )
add_subdirectory( AssemblerLib )
add_subdirectory( BaseLib )
# TODO This is a hack but we have to make sure that Boost is built first
if(TARGET Boost)
    add_dependencies(BaseLib Boost)
endif()
add_subdirectory( DataHolderLib )
add_subdirectory( FileIO )
add_subdirectory( GeoLib )
add_subdirectory( InSituLib )
add_subdirectory( MaterialsLib/Adsorption )
add_subdirectory( MathLib )
add_subdirectory( MeshLib )
add_subdirectory( MeshGeoToolsLib )
add_subdirectory( NumLib )
add_subdirectory( ProcessLib )
if( OGS_BUILD_TESTS AND NOT IS_SUBPROJECT )
    add_subdirectory( Tests )

    if(OGS_USE_MPI)
        add_subdirectory( SimpleTests/MeshTests/MPI )
    else()
        add_subdirectory( SimpleTests/MatrixTests )
        add_subdirectory( SimpleTests/MeshTests )
        if(NOT MSVC AND BLAS_FOUND AND LAPACK_FOUND)
            add_subdirectory( SimpleTests/SolverTests )
        endif()
    endif()
endif() # OGS_BUILD_TESTS

# The configuration must be called from the source dir and not BaseLib/.
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/BaseLib/BuildInfo.cpp.in"
    "${CMAKE_CURRENT_BINARY_DIR}/BaseLib/BuildInfo.cpp" @ONLY)

include(scripts/cmake/MarkVariablesAdvanced.cmake)
