# LIE; ThermoThermoHydroMechanics

AddTest(
    NAME LARGE_LIE_THM_TaskB_HM
    PATH LIE/ThermoHydroMechanics
    RUNTIME 60
    EXECUTABLE ogs
    EXECUTABLE_ARGS TaskB_HM.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    DIFF_DATA
    expected_TaskB_HM_pcs_0_ts_4_t_18.000000.vtu TaskB_HM_pcs_0_ts_4_t_18.000000.vtu pressure pressure 1e-12 1e-12
    expected_TaskB_HM_pcs_0_ts_4_t_18.000000.vtu TaskB_HM_pcs_0_ts_4_t_18.000000.vtu pressure_interpolated pressure_interpolated 1e-12 1e-12
    expected_TaskB_HM_pcs_0_ts_4_t_18.000000.vtu TaskB_HM_pcs_0_ts_4_t_18.000000.vtu displacement displacement 1e-12 1e-12
    expected_TaskB_HM_pcs_0_ts_4_t_18.000000.vtu TaskB_HM_pcs_0_ts_4_t_18.000000.vtu displacement_jump1 displacement_jump1 1e-12 1e-12
    expected_TaskB_HM_pcs_0_ts_4_t_18.000000.vtu TaskB_HM_pcs_0_ts_4_t_18.000000.vtu nodal_w nodal_w 1e-12 1e-12
    expected_TaskB_HM_pcs_0_ts_4_t_18.000000.vtu TaskB_HM_pcs_0_ts_4_t_18.000000.vtu nodal_aperture nodal_aperture 1e-12 1e-12
    expected_TaskB_HM_pcs_0_ts_4_t_18.000000.vtu TaskB_HM_pcs_0_ts_4_t_18.000000.vtu strain_xx strain_xx 1e-12 1e-12
    expected_TaskB_HM_pcs_0_ts_4_t_18.000000.vtu TaskB_HM_pcs_0_ts_4_t_18.000000.vtu strain_yy strain_yy 1e-12 1e-12
    expected_TaskB_HM_pcs_0_ts_4_t_18.000000.vtu TaskB_HM_pcs_0_ts_4_t_18.000000.vtu strain_xy strain_xy 1e-12 1e-12
    expected_TaskB_HM_pcs_0_ts_4_t_18.000000.vtu TaskB_HM_pcs_0_ts_4_t_18.000000.vtu stress_xx stress_xx 1e-4 1e-10
    expected_TaskB_HM_pcs_0_ts_4_t_18.000000.vtu TaskB_HM_pcs_0_ts_4_t_18.000000.vtu stress_yy stress_yy 1e-4 1e-10
    expected_TaskB_HM_pcs_0_ts_4_t_18.000000.vtu TaskB_HM_pcs_0_ts_4_t_18.000000.vtu velocity velocity 1e-12 1e-12
)

AddTest(
    NAME LARGE_LIE_THM_Lauwerier
    PATH LIE/ThermoHydroMechanics
    RUNTIME 60
    EXECUTABLE ogs
    EXECUTABLE_ARGS Lauwerier.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    DIFF_DATA
    expected_Lauwerier_pcs_0_ts_2_t_4000000.000000.vtu Lauwerier_pcs_0_ts_2_t_4000000.000000.vtu pressure_interpolated pressure_interpolated 1e-12 1e-12
    expected_Lauwerier_pcs_0_ts_2_t_4000000.000000.vtu Lauwerier_pcs_0_ts_2_t_4000000.000000.vtu temperature_interpolated temperature_interpolated 1e-12 1e-12
    expected_Lauwerier_pcs_0_ts_2_t_4000000.000000.vtu Lauwerier_pcs_0_ts_2_t_4000000.000000.vtu velocity velocity 1e-12 1e-12
)

AddTest(
    NAME LIE_THM_single_crack_inside_nomatflow
    PATH LIE/ThermoHydroMechanics
    RUNTIME 60
    EXECUTABLE ogs
    EXECUTABLE_ARGS single_crack_inside_nomatflow.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    DIFF_DATA
    expected_single_crack_inside_nomatflow_pcs_0_ts_10_t_1.000000.vtu single_crack_inside_nomatflow_pcs_0_ts_10_t_1.000000.vtu pressure_interpolated pressure_interpolated 1e-12 1e-12
    expected_single_crack_inside_nomatflow_pcs_0_ts_10_t_1.000000.vtu single_crack_inside_nomatflow_pcs_0_ts_10_t_1.000000.vtu temperature_interpolated temperature_interpolated 1e-12 1e-12
    expected_single_crack_inside_nomatflow_pcs_0_ts_10_t_1.000000.vtu single_crack_inside_nomatflow_pcs_0_ts_10_t_1.000000.vtu displacement_jump1 displacement_jump1 1e-12 1e-12
)


AddTest(
    NAME LIE_THM_single_crack_inside_matflow
    PATH LIE/ThermoHydroMechanics
    RUNTIME 60
    EXECUTABLE ogs
    EXECUTABLE_ARGS single_crack_inside_matflow.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    DIFF_DATA
    expected_single_crack_inside_matflow_pcs_0_ts_10_t_1.000000.vtu single_crack_inside_matflow_pcs_0_ts_10_t_1.000000.vtu pressure_interpolated pressure_interpolated 1e-12 1e-12
    expected_single_crack_inside_matflow_pcs_0_ts_10_t_1.000000.vtu single_crack_inside_matflow_pcs_0_ts_10_t_1.000000.vtu temperature_interpolated temperature_interpolated 1e-12 1e-12
    expected_single_crack_inside_matflow_pcs_0_ts_10_t_1.000000.vtu single_crack_inside_matflow_pcs_0_ts_10_t_1.000000.vtu displacement_jump1 displacement_jump1 1e-12 1e-12
)


AddTest(
    NAME LIE_THM_two_cracks_branch
    PATH LIE/ThermoHydroMechanics
    RUNTIME 60
    EXECUTABLE ogs
    EXECUTABLE_ARGS two_cracks_branch.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    DIFF_DATA
    expected_two_cracks_branch_pcs_0_ts_10_t_1.000000.vtu two_cracks_branch_pcs_0_ts_10_t_1.000000.vtu pressure_interpolated pressure_interpolated 1e-12 1e-12
    expected_two_cracks_branch_pcs_0_ts_10_t_1.000000.vtu two_cracks_branch_pcs_0_ts_10_t_1.000000.vtu temperature_interpolated temperature_interpolated 1e-12 1e-12
    expected_two_cracks_branch_pcs_0_ts_10_t_1.000000.vtu two_cracks_branch_pcs_0_ts_10_t_1.000000.vtu displacement_jump1 displacement_jump1 1e-12 1e-12
    expected_two_cracks_branch_pcs_0_ts_10_t_1.000000.vtu two_cracks_branch_pcs_0_ts_10_t_1.000000.vtu displacement_jump2 displacement_jump2 1e-12 1e-12
)

AddTest(
    NAME LIE_THM_two_cracks_junction
    PATH LIE/ThermoHydroMechanics
    RUNTIME 60
    EXECUTABLE ogs
    EXECUTABLE_ARGS two_cracks_junction.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
    DIFF_DATA
    expected_two_cracks_junction_pcs_0_ts_10_t_1.000000.vtu two_cracks_junction_pcs_0_ts_10_t_1.000000.vtu pressure_interpolated pressure_interpolated 1e-12 1e-12
    expected_two_cracks_junction_pcs_0_ts_10_t_1.000000.vtu two_cracks_junction_pcs_0_ts_10_t_1.000000.vtu temperature_interpolated temperature_interpolated 1e-12 1e-12
    expected_two_cracks_junction_pcs_0_ts_10_t_1.000000.vtu two_cracks_junction_pcs_0_ts_10_t_1.000000.vtu displacement_jump1 displacement_jump1 1e-12 1e-12
    expected_two_cracks_junction_pcs_0_ts_10_t_1.000000.vtu two_cracks_junction_pcs_0_ts_10_t_1.000000.vtu displacement_jump2 displacement_jump2 1e-12 1e-12
    expected_two_cracks_junction_pcs_0_ts_10_t_1.000000.vtu two_cracks_junction_pcs_0_ts_10_t_1.000000.vtu displacement_jump3 displacement_jump3 1e-12 1e-12
)
