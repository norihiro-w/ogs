#Source files grouped by a directory
GET_SOURCE_FILES(SOURCES_SOLUTIONLIB)
set(SOURCES ${SOURCES_SOLUTIONLIB})
GET_SOURCE_FILES(SOURCES_CORE Core)
set(SOURCES ${SOURCES} ${SOURCES_CORE})
GET_SOURCE_FILES(SOURCES_FEM Fem)
set(SOURCES ${SOURCES} ${SOURCES_FEM})


include_directories(
	.
	${CMAKE_CURRENT_SOURCE_DIR}/../BaseLib
	${CMAKE_CURRENT_SOURCE_DIR}/../GeoLib
	${CMAKE_CURRENT_SOURCE_DIR}/../MathLib
	${CMAKE_CURRENT_SOURCE_DIR}/../MeshLib
	${CMAKE_CURRENT_SOURCE_DIR}/../MeshGeoToolsLib
	${CMAKE_CURRENT_SOURCE_DIR}/../NumLib
)

# Create the library
add_library(SolutionLib STATIC ${SOURCES})

set_target_properties(SolutionLib PROPERTIES LINKER_LANGUAGE CXX)

target_link_libraries(SolutionLib
	BaseLib
	GeoLib
	MathLib
	MeshLib
	MeshGeoToolsLib
	AssemblerLib
	NumLib
	logog
)

