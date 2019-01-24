# AddTest(
#     NAME SmallDeformationWithPT_3D_ThermoElastic_Stress_Analysis
#     PATH SmallDeformationWithPT
#     EXECUTABLE ogs
#     EXECUTABLE_ARGS cube_1e3.prj
#     WRAPPER time
#     TESTER vtkdiff
#     REQUIREMENTS NOT (OGS_USE_LIS OR OGS_USE_MPI)
#     RUNTIME 17
#     DIFF_DATA
#     stress_analytical.vtu cube_1e3_tm_pcs_0_ts_17_t_72000.000000.vtu sigma sigma 1e-5 1e-12
#     expected_cube_1e3_tm_pcs_0_ts_17_t_72000.000000.vtu cube_1e3_tm_pcs_0_ts_17_t_72000.000000.vtu displacement displacement 1e-10 1e-12
#     expected_cube_1e3_tm_pcs_0_ts_17_t_72000.000000.vtu cube_1e3_tm_pcs_0_ts_17_t_72000.000000.vtu temperature temperature 1e-10 1e-12
#     expected_cube_1e3_tm_pcs_0_ts_17_t_72000.000000.vtu cube_1e3_tm_pcs_0_ts_17_t_72000.000000.vtu sigma sigma 1e-6 1e-12
#     expected_cube_1e3_tm_pcs_0_ts_17_t_72000.000000.vtu cube_1e3_tm_pcs_0_ts_17_t_72000.000000.vtu epsilon epsilon 1e-16 0
# )

