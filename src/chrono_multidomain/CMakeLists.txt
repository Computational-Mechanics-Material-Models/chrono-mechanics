#=============================================================================
# CMake configuration file for the Chrono multidomain module
#=============================================================================

option(CH_ENABLE_MODULE_MULTIDOMAIN "Enable the Chrono multidomain module" OFF)

if(NOT CH_ENABLE_MODULE_MULTIDOMAIN)
    return()
endif()

message(STATUS "\n==== Chrono multidomain module ====\n")

if(NOT MPI_CXX_FOUND)
    message("  Chrono::Multidomain requires MPI for C++, but MPI was not found. Disabling Chrono::Multidomain.")
    set(CH_ENABLE_MODULE_MULTIDOMAIN OFF CACHE BOOL "Enable the Chrono multidomain module" FORCE)
    return()
endif()

set(Chrono_MULTIDOMAIN_SOURCES
    ChDomain.cpp
    ChDomainBuilder.cpp
    ChDomainManager.cpp
    ChDomainManagerSharedmemory.cpp
    ChDomainManagerMPI.cpp
    ChSystemDescriptorMultidomain.cpp
    ChSolverPSORmultidomain.cpp
    ChSolverLumpedMultidomain.cpp
    ChSolverBBmultidomain.cpp
    ChMpi.cpp
)

set(Chrono_MULTIDOMAIN_HEADERS
    ChApiMultiDomain.h
    ChDomain.h
    ChDomainBuilder.h
    ChDomainManager.h
    ChDomainManagerSharedmemory.h
    ChDomainManagerMPI.h
    ChSystemDescriptorMultidomain.h
    ChSolverPSORmultidomain.h
    ChSolverLumpedMultidomain.h
    ChSolverBBmultidomain.h
    ChMpi.h
)

source_group("" FILES ${Chrono_MULTIDOMAIN_SOURCES} ${Chrono_MULTIDOMAIN_HEADERS})

add_library(Chrono_multidomain
            ${Chrono_MULTIDOMAIN_SOURCES}
            ${Chrono_MULTIDOMAIN_HEADERS})
add_library(Chrono::multidomain ALIAS Chrono_multidomain)

if(CH_WHOLE_PROG_OPT)
  set_target_properties(Chrono_multidomain PROPERTIES COMPILE_FLAGS "/GL")
  set_target_properties(Chrono_multidomain PROPERTIES LINK_FLAGS "/LTCG")
endif()

if(CH_STATIC)
  set_target_properties(Chrono_multidomain PROPERTIES POSITION_INDEPENDENT_CODE ON)
endif()

if(MSVC)
  set_target_properties(Chrono_multidomain PROPERTIES MSVC_RUNTIME_LIBRARY ${CH_MSVC_RUNTIME_LIBRARY})
endif()

set_target_properties(Chrono_multidomain PROPERTIES
                      COMPILE_FLAGS "${CH_CXX_FLAGS} ${MPI_CXX_COMPILE_FLAGS}"
                      LINK_FLAGS "${CH_LINKERFLAG_LIB} ${MPI_CXX_LINK_FLAGS}")

target_include_directories(Chrono_multidomain SYSTEM PRIVATE ${MPI_CXX_INCLUDE_PATH})

target_compile_definitions(Chrono_multidomain PRIVATE "CH_API_COMPILE_MULTIDOMAIN")
target_compile_definitions(Chrono_multidomain PRIVATE "CH_IGNORE_DEPRECATED")

target_link_libraries(Chrono_multidomain
                      Chrono_core
                      ${MPI_CXX_LIBRARIES})

install(TARGETS Chrono_multidomain
  EXPORT ChronoTargets
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  INCLUDES DESTINATION include/chrono_multidomain)

install(FILES ${Chrono_MULTIDOMAIN_HEADERS}
        DESTINATION include/chrono_multidomain)

set(DEPENDENCIES_MULTIDOMAIN "" PARENT_SCOPE)
