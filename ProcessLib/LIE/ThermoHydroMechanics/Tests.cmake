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
