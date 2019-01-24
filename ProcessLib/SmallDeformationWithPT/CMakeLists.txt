APPEND_SOURCE_FILES(SOURCES)

add_library(SmallDeformationWithPT ${SOURCES})
if(BUILD_SHARED_LIBS)
    install(TARGETS SmallDeformationWithPT LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()
target_link_libraries(SmallDeformationWithPT PUBLIC ProcessLib)

include(Tests.cmake)
