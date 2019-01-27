APPEND_SOURCE_FILES(SOURCES)

add_library(THM ${SOURCES})
if(BUILD_SHARED_LIBS)
    install(TARGETS THM LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()
target_link_libraries(THM PUBLIC ProcessLib)

include(Tests.cmake)
