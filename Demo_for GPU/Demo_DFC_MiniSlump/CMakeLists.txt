cmake_minimum_required(VERSION 3.16.3)
cmake_policy(SET CMP0091 NEW)

project(MiniSlump LANGUAGES CXX)

# The local Chrono 10 build contains Release libraries only.
if(CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_CONFIGURATION_TYPES Release CACHE STRING "Available build configurations" FORCE)
elseif(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

find_package(Chrono
  REQUIRED
  COMPONENTS Irrlicht PardisoMKL
  CONFIG
)

set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(CHRONO_CXX_STANDARD MATCHES "^cxx_std_([0-9]+)$")
  set(CMAKE_CXX_STANDARD "${CMAKE_MATCH_1}")
else()
  set(CMAKE_CXX_STANDARD "${CHRONO_CXX_STANDARD}")
endif()

if(APPLE)
  set(CMAKE_MACOSX_BUNDLE ON)
endif()

add_executable(MiniSlump
  demo_MiniSlump.cpp
  FreshConcreteContact_V10.cpp
)

# FreshConcreteContact_V10.cpp is included directly by demo_MiniSlump.cpp.
# Listing it in the project is useful for IDEs, but it must not be compiled
# as a separate translation unit.
set_source_files_properties(FreshConcreteContact_V10.cpp PROPERTIES HEADER_FILE_ONLY TRUE)

target_include_directories(MiniSlump PRIVATE ${CHRONO_INCLUDE_DIRS})

if(DEFINED ENV{CONDA_PREFIX} AND EXISTS "$ENV{CONDA_PREFIX}/Library/include/boost/numeric/odeint.hpp")
  target_include_directories(MiniSlump PRIVATE "$ENV{CONDA_PREFIX}/Library/include")
elseif(EXISTS "C:/Users/Jiaqi/anaconda3/Library/include/boost/numeric/odeint.hpp")
  target_include_directories(MiniSlump PRIVATE "C:/Users/Jiaqi/anaconda3/Library/include")
endif()

target_compile_definitions(MiniSlump PRIVATE
  "CHRONO_DATA_DIR=\"${CHRONO_DATA_DIR}\""
)
target_compile_options(MiniSlump PRIVATE ${CHRONO_CXX_FLAGS})
target_link_options(MiniSlump PRIVATE ${CHRONO_LINKER_FLAGS})

if(MSVC)
  target_compile_definitions(MiniSlump PRIVATE
    _ENABLE_EXTENDED_ALIGNED_STORAGE
    _CRT_SECURE_NO_DEPRECATE
    _SCL_SECURE_NO_DEPRECATE
    NOMINMAX
  )
endif()

if(CHRONO_TARGETS)
  target_link_libraries(MiniSlump PRIVATE ${CHRONO_TARGETS})
else()
  target_link_libraries(MiniSlump PRIVATE ${CHRONO_LIBRARIES})
endif()

# Copy simulation inputs beside the executable. The application resolves these
# paths relative to its own location.
foreach(_input_file
    minislump-Cut4.obj
    DFCgeo000-data-particles-mini-2x4-dist.dat
)
  add_custom_command(TARGET MiniSlump POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
            "${CMAKE_SOURCE_DIR}/${_input_file}"
            "$<TARGET_FILE_DIR:MiniSlump>/${_input_file}"
  )
endforeach()

if(COMMAND add_DLL_copy_command)
  add_DLL_copy_command()
endif()

# Chrono 10 package configurations do not always provide the DLL copy helper.
if(DEFINED Chrono_DIR)
  get_filename_component(_chrono_cmake_dir "${Chrono_DIR}" ABSOLUTE)
  get_filename_component(_chrono_build_dir "${_chrono_cmake_dir}/.." ABSOLUTE)
endif()

set(_extra_runtime_dlls
  "${_chrono_build_dir}/bin/Release/Chrono_core.dll"
  "${_chrono_build_dir}/bin/Release/Chrono_irrlicht.dll"
  "${_chrono_build_dir}/bin/Release/Chrono_pardisomkl.dll"
  "C:/Users/Jiaqi/library/library/irrlicht-1.8.5/bin/Win64-VisualStudio/Irrlicht.dll"
  "C:/Windows/System32/VCOMP140.DLL"
)

foreach(_runtime_dll IN LISTS _extra_runtime_dlls)
  if(EXISTS "${_runtime_dll}")
    add_custom_command(TARGET MiniSlump POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy_if_different
              "${_runtime_dll}"
              "$<TARGET_FILE_DIR:MiniSlump>"
    )
  endif()
endforeach()

# PardisoMKL depends on the Intel MKL and OpenMP redistributables.
file(GLOB _intel_runtime_dlls
  "C:/Program Files (x86)/Intel/oneAPI/mkl/latest/bin/mkl*.2.dll"
)
list(APPEND _intel_runtime_dlls
  "C:/Program Files (x86)/Intel/oneAPI/compiler/latest/bin/libiomp5md.dll"
)

foreach(_intel_dll IN LISTS _intel_runtime_dlls)
  if(EXISTS "${_intel_dll}")
    add_custom_command(TARGET MiniSlump POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy_if_different
              "${_intel_dll}"
              "$<TARGET_FILE_DIR:MiniSlump>"
    )
  endif()
endforeach()
