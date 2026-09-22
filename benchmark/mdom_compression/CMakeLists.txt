cmake_minimum_required(VERSION 3.18)

project(chrono_mdom_compression LANGUAGES CXX)

find_package(MPI REQUIRED COMPONENTS CXX)
find_package(Chrono REQUIRED COMPONENTS LDPM MULTIDOMAIN CONFIG)

add_executable(ldpm_uc LDPM_UC.cpp)
target_compile_features(ldpm_uc PRIVATE cxx_std_17)
target_compile_definitions(
  ldpm_uc
  PRIVATE
    "CHRONO_DATA_DIR=\"${CHRONO_DATA_DIR}\""
    "LDPM_EXAMPLE_DIR=\"${CMAKE_CURRENT_SOURCE_DIR}/../../template_project_ldpm/Unconfined Compression\"")
target_link_libraries(
  ldpm_uc
  PRIVATE
    Chrono::Chrono_core
    Chrono::Chrono_ldpm
    Chrono::Chrono_multidomain
    MPI::MPI_CXX)
