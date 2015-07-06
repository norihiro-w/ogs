#-------------------------------------------------------------------------------
# CMake script to configure a module(s) in a directory
#
# If this directory contains sub directories including modules, please write
#
#  SET (MODULE_SUBDIRECTORIES TRUE)
#
# and skip the followings.
# 
# To set up a single module, this script should set the following variables
# - MODULE_PROCESS_HEADERS
#     Path to header files defining processes
# - MODULE_PROCESS_LIST
#     Processs registration list. Please use the following commands to register
#     processes. 
#     - OGS_ADD_PROCESS
#     - OGS_ADD_PROCESS_SYS
#     - OGS_ADD_PROCESS_SYS_SOLVER
# - MODULE_TEST_INPUT_DIR (optional)
#     Path to directories where test input data are stored
# - MODULE_SOURCE_FILES (optional)
# - MODULE_TEST_FILES (optional)
# 
# The following variables are available to configure the above.
# - MODULE_DIR  : a path to this module
# - MODULE_NAME : This module name, i.e. directory name 
#-------------------------------------------------------------------------------

# Set process header files
SET (MODULE_PROCESS_HEADERS
    FunctionTemperature.h
)

# Register processes
OGS_ADD_PROCESS( MODULE_PROCESS_LIST 
    "HEAT_TRANSPORT" 
    FunctionTemperature 
    )

# Set test directories
SET (MODULE_TEST_INPUT_DIR
    ogata-banks
)
