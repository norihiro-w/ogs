include(${PROJECT_SOURCE_DIR}/scripts/cmake/OGSEnabledElements.cmake)

add_subdirectory(ApplicationsLib)
add_subdirectory(DataHolderLib)
add_subdirectory(FileIO)

add_subdirectory(Utils/PostProcessing)
AddTest(
    NAME postLIE
    PATH LIE/PostProcessing
    EXECUTABLE postLIE
    EXECUTABLE_ARGS -i single_joint_pcs_0.pvd -o ${Data_BINARY_DIR}/LIE/PostProcessing/post_single_joint_pcs_0.pvd
    REQUIREMENTS NOT OGS_USE_MPI AND OGS_BUILD_PROCESS_LIE
    TESTER vtkdiff
    DIFF_DATA
    expected_post_single_joint_pcs_0_ts_1_t_1.000000.vtu post_single_joint_pcs_0_ts_1_t_1.000000.vtu u u 1e-14 1e-14
)

if(OGS_BUILD_UTILS AND NOT IS_SUBPROJECT)
    add_subdirectory(Utils)
elseif(OGS_BUILD_GUI)
    add_subdirectory(Utils/OGSFileConverter)
endif() # OGS_BUILD_UTILS AND NOT IS_SUBPROJECT

if(OGS_BUILD_GUI)
    add_subdirectory(DataExplorer)
endif() # OGS_BUILD_GUI

if(OGS_BUILD_CLI)
    add_subdirectory(CLI)
endif() # OGS_BUILD_CLI

if(OGS_INSITU)
    add_subdirectory(InSituLib)
endif()
