#===============================================================================
# CMake configuration file for the Chrono dfc module
#===============================================================================

option(ENABLE_MODULE_DFC "Enable the Chrono DFC module" OFF)

# Return now if this module is not enabled
IF(NOT ENABLE_MODULE_DFC)
  RETURN()
ENDIF()

MESSAGE(STATUS "\n==== Chrono DFC module ====\n")

# List all source and header files for the ChronoEngine_DFC library
SET(ChronoEngine_DFC_SOURCES 
    FreshConcreteContact.cpp
)

SET(ChronoEngine_DFC_HEADERS

)

SOURCE_GROUP("" FILES 
            ${ChronoEngine_DFC_SOURCES} 
            ${ChronoEngine_DFC_HEADERS}
)

# Include directories (ensure DFC_INCLUDE_DIR is defined or specify the correct path)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

# Create the ChronoEngine_DFC library
ADD_LIBRARY(ChronoEngine_dfc
  ${ChronoEngine_DFC_SOURCES}
  ${ChronoEngine_DFC_HEADERS}
)

# Set compiler and linker flags
SET_TARGET_PROPERTIES(ChronoEngine_dfc PROPERTIES
  COMPILE_FLAGS "${CH_CXX_FLAGS}"
  LINK_FLAGS "${CH_LINKERFLAG_LIB}"
)

# Compile definitions
target_compile_definitions(ChronoEngine_dfc PRIVATE "CH_DFC_API_COMPILE")
target_compile_definitions(ChronoEngine_dfc PRIVATE "CH_IGNORE_DEPRECATED")

# Link required libraries
TARGET_LINK_LIBRARIES(ChronoEngine_dfc
  ChronoEngine
  ${DFC_LIBRARY}
)


# Add dependencies
ADD_DEPENDENCIES(ChronoEngine_dfc ChronoEngine)

# Handle Windows DLL copying
#IF(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
#  ADD_CUSTOM_COMMAND(
#    TARGET ChronoEngine_dfc POST_BUILD
#    COMMAND ${CMAKE_COMMAND} -E make_directory
#            "${CMAKE_BINARY_DIR}/bin/$<CONFIGURATION>"
#    COMMAND ${CMAKE_COMMAND} -E copy_if_different
#            "${CH_dfc_DLL}"
#            "${CMAKE_BINARY_DIR}/bin/$<CONFIGURATION>"
#    COMMENT "Copying dfc DLL..."
#  )
#  install(FILES "${CH_dfc_DLL}" DESTINATION bin)
#ENDIF()

# Install the library and headers
INSTALL(TARGETS ChronoEngine_dfc
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)

INSTALL(FILES ${ChronoEngine_DFC_HEADERS} DESTINATION include/chrono_dfc)

