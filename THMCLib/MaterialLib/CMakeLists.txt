#Source files grouped by a directory
# ./
FILE(GLOB MaterialLib_Files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} MaterialLib/*.h MaterialLib/*.cpp)
SOURCE_GROUP( ogs6THMC\\MaterialLib FILES ${MaterialLib_Files} )
# ./Core
#FILE(GLOB DiscreteLib_Core_Files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} DiscreteLib/Core/*.h DiscreteLib/Core/*.cpp)
#SOURCE_GROUP( ogsNumerics\\DiscreteLib\\Core FILES ${DiscreteLib_Core_Files} )


#
#LIST(APPEND MaterialLib_Files
#    ${MaterialLib_Files} 
#) 

