
include_directories(
	${CMAKE_SOURCE_DIR}/Utils/FileConverter
	${CMAKE_SOURCE_DIR}/BaseLib
	${CMAKE_SOURCE_DIR}/GeoLib
	${CMAKE_SOURCE_DIR}/FileIO
	${CMAKE_SOURCE_DIR}/MathLib
	${CMAKE_SOURCE_DIR}/MeshLib
	${CMAKE_SOURCE_DIR}/Gui/OpenSG
)

if (QT4_FOUND)
	add_executable(ConvertSHPToGLI ConvertSHPToGLI.cpp)
	set_target_properties(ConvertSHPToGLI PROPERTIES FOLDER Utilities)

	target_link_libraries(ConvertSHPToGLI
		FileIO
		GeoLib
		BaseLib
		shp
		${QT_LIBRARIES}
	)

	add_executable(generateBCandGLI generateBCandGLI.cpp )
	target_link_libraries(generateBCandGLI
		ApplicationsLib
		FileIO
		GeoLib
		${QT_LIBRARIES}
	)

	add_executable(generateBCFromPolyline generateBCFromPolyline.cpp )
	target_link_libraries(generateBCFromPolyline
		ApplicationsLib
		FileIO
		GeoLib
		${QT_LIBRARIES}
	)

	set_target_properties(generateBCandGLI generateBCFromPolyline
		PROPERTIES FOLDER Utilities)

	add_executable(FEFLOW2OGS FEFLOW2OGS.cpp)
	set_target_properties(FEFLOW2OGS PROPERTIES FOLDER Utilities)
	target_link_libraries(FEFLOW2OGS
		GeoLib
		MeshLib
		FileIO
		InSituLib
		${CATALYST_LIBRARIES}
		${QT_LIBRARIES}
	)
	ADD_CATALYST_DEPENDENCY(FEFLOW2OGS)
endif () # QT4_FOUND

if (Boost_FOUND)
	add_executable(generateMatPropsFromMatID generateMatPropsFromMatID.cpp )
	target_link_libraries(generateMatPropsFromMatID
		FileIO
		MeshLib
		logog
		InSituLib
		${CATALYST_LIBRARIES}
	)
	ADD_CATALYST_DEPENDENCY(generateMatPropsFromMatID)
	set_target_properties(generateMatPropsFromMatID
		PROPERTIES FOLDER Utilities)
endif () # Boost_FOUND

add_executable(GMSH2OGS GMSH2OGS.cpp)
set_target_properties(GMSH2OGS PROPERTIES FOLDER Utilities)
target_link_libraries(GMSH2OGS
	MeshLib
	FileIO
	InSituLib
	${CATALYST_LIBRARIES}
)
ADD_CATALYST_DEPENDENCY(GMSH2OGS)
if(MSVC)
	target_link_libraries(GMSH2OGS winmm)
endif()

add_executable(OGS2VTK OGS2VTK.cpp)
set_target_properties(OGS2VTK PROPERTIES FOLDER Utilities)
target_link_libraries(OGS2VTK
	MeshLib
	FileIO
	zlib
	InSituLib
	${CATALYST_LIBRARIES}
)
ADD_CATALYST_DEPENDENCY(OGS2VTK)

add_executable(VTK2OGS VTK2OGS.cpp)
set_target_properties(VTK2OGS PROPERTIES FOLDER Utilities)
target_link_libraries(VTK2OGS
	MeshLib
	FileIO
	zlib
	InSituLib
	${CATALYST_LIBRARIES}
)
ADD_CATALYST_DEPENDENCY(VTK2OGS)

add_executable(VTK2TIN VTK2TIN.cpp)
set_target_properties(VTK2TIN PROPERTIES FOLDER Utilities)
target_link_libraries(VTK2TIN
	MeshLib
	FileIO
	zlib
	InSituLib
	${CATALYST_LIBRARIES}
)
ADD_CATALYST_DEPENDENCY(VTK2TIN)

add_executable(TIN2VTK TIN2VTK.cpp)
set_target_properties(TIN2VTK PROPERTIES FOLDER Utilities)
target_link_libraries(TIN2VTK
	MeshLib
	FileIO
	zlib
	InSituLib
	${CATALYST_LIBRARIES}
)
ADD_CATALYST_DEPENDENCY(TIN2VTK)

####################
### Installation ###
####################
install(TARGETS generateMatPropsFromMatID GMSH2OGS OGS2VTK VTK2OGS VTK2TIN
	RUNTIME DESTINATION bin COMPONENT ogs_converter)

if(QT4_FOUND)
	install(TARGETS ConvertSHPToGLI generateBCandGLI generateBCFromPolyline generateBCandGLI
		FEFLOW2OGS RUNTIME DESTINATION bin COMPONENT ogs_converter)
endif()

cpack_add_component(ogs_converter
	DISPLAY_NAME "File converter"
	DESCRIPTION "File converter."
	GROUP Utilities
)

ADD_EXECUTABLE (Raster2Mesh Raster2Mesh.cpp)
SET_TARGET_PROPERTIES(Raster2Mesh PROPERTIES FOLDER Utilities)
TARGET_LINK_LIBRARIES (Raster2Mesh
	MeshLib
	FileIO
	zlib
)

ADD_EXECUTABLE (Raster2QuadMesh Raster2QuadMesh.cpp)
SET_TARGET_PROPERTIES(Raster2QuadMesh PROPERTIES FOLDER Utilities)
TARGET_LINK_LIBRARIES (Raster2QuadMesh
	MeshLib
	FileIO
	zlib
)

ADD_EXECUTABLE (CSV2Raster CSV2Raster.cpp)
SET_TARGET_PROPERTIES(CSV2Raster PROPERTIES FOLDER Utilities)
TARGET_LINK_LIBRARIES (CSV2Raster
	MeshLib
	FileIO
	zlib
)

ADD_EXECUTABLE (addMarginToRaster addMarginToRaster.cpp)
SET_TARGET_PROPERTIES(addMarginToRaster PROPERTIES FOLDER Utilities)
TARGET_LINK_LIBRARIES (addMarginToRaster
	MeshLib
	FileIO
	zlib
)
