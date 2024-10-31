# Set minimum version of CMake.
cmake_minimum_required(VERSION 3.12)

# Set project name and version
project(dhb VERSION 1.0.0 LANGUAGES CXX)

include(GNUInstallDirs)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


# Configure the option whether you want to use the systems own allocater
# DHB_WITH_SYSTEM_ALLOCATOR=ON
# or our own written allocator strategy
# DHB_WITH_SYSTEM_ALLOCATOR=OFF
option(DHB_WITH_SYSTEM_ALLOCATOR "Use DHB System Allocator" OFF)
option(DHB_WITH_64BIT_IDS "Use 64 bit IDs." OFF)
option(DHB_WARNINGS_AS_ERRORS "Sets warning level high, treats warnings as errors ." OFF)

add_library(dhb STATIC)

set(public_headers
  include/dhb/vec.h
  include/dhb/batcher.h
  include/dhb/block.h
  include/dhb/buckets.h
  include/dhb/dynamic_hashed_blocks.h
  include/dhb/graph.h
  include/dhb/hash_tools.h
  include/dhb/integer_log2.h
  include/dhb/submatrix.h
)

target_sources(${PROJECT_NAME}
  PRIVATE
    src/buckets.cpp
    src/dynamic_hashed_blocks.cpp
    src/graph.cpp
)

# ===================================
# Options check and set definitions
# ===================================
if (DHB_WITH_SYSTEM_ALLOCATOR)
  MESSAGE(STATUS "DHB uses system allocator.")
  target_compile_definitions(dhb PUBLIC -DDHB_SYSTEM_ALLOCATOR)
endif()

if (DHB_WITH_64BIT_IDS)
  MESSAGE(STATUS "DHB uses 64 bit vertex IDs.")
  target_compile_definitions(dhb PUBLIC -DDHB_64BIT_IDS)
else()
  MESSAGE(STATUS "DHB uses 32 bit vertex IDs.")
endif()

if (DHB_WARNINGS_AS_ERRORS)
  MESSAGE(STATUS "DHB is extra pedantic, throws errors on warnings.")
  target_compile_options(dhb PRIVATE -Werror -Wall -Wextra -Wpedantic -Wmaybe-uninitialized)
endif()

find_package(OpenMP)
# Thanks to NetworKit for the following OpenMP code looking or OpenMP in the
# most common folders:
# https://github.com/networkit/networkit/blob/master/CMakeLists.txt#L145
if(NOT OpenMP_FOUND)
  if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang")
      # This will find default libomp-installations for homebrew/MacPorts
      find_library(LIBOMP_PATH NAMES omp PATHS "/usr/local/opt/libomp/lib" "/opt/local/lib/libomp" "/opt/homebrew/opt/libomp/lib")
      find_path(LIBOMP_INCLUDE NAMES omp.h PATHS "/usr/local/opt/libomp/include" "/opt/local/include/libomp" "/opt/homebrew/opt/libomp/include")
      if(LIBOMP_PATH AND LIBOMP_INCLUDE)
          set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp -I${LIBOMP_INCLUDE}" CACHE STRING "Manually set" FORCE)
      endif()
  elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
      if(DEFINED ENV{CONDA_PREFIX})
          find_library(LIBOMP_PATH NAMES omp HINTS "$ENV{CONDA_PREFIX}/lib"
              PATHS "/usr/lib" "/usr/lib64")
          find_path(LIBOMP_INCLUDE NAMES omp.h HINTS "$ENV{CONDA_PREFIX}/include"
              PATHS "/usr/include")
      else()
          find_library(LIBOMP_PATH NAMES omp PATHS "/usr/lib" "/usr/lib64")
          find_path(LIBOMP_INCLUDE NAMES omp.h PATHS "/usr/include")
      endif()
      if(LIBOMP_PATH AND LIBOMP_INCLUDE)
          set(OpenMP_CXX_FLAGS "-fopenmp -I${LIBOMP_INCLUDE}" CACHE STRING "Manually set" FORCE)
      endif()
  endif()

  # Set OpenMP-folders in case they are found with the aid of additional hints
  if(LIBOMP_PATH AND LIBOMP_INCLUDE)
      set(OpenMP_CXX_LIB_NAMES "omp" CACHE STRING "Manually set" FORCE)
      set(OpenMP_omp_LIBRARY "${LIBOMP_PATH}" CACHE STRING "Manually set" FORCE)
  else()
      message(FATAL_ERROR "libomp was not found, but necessary to configure DHB with ${CMAKE_CXX_COMPILER_ID}")
  endif()

  # After setting basic OpenMP-folders, run find_package again to set everything. Also acts as a final sanity check.
  find_package(OpenMP REQUIRED)
endif()
target_link_libraries(${PROJECT_NAME} PUBLIC OpenMP::OpenMP_CXX)

target_include_directories(${PROJECT_NAME}
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

install(TARGETS ${PROJECT_NAME}
    EXPORT "${PROJECT_NAME}Targets"
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(FILES ${public_headers} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME})

install(EXPORT "${PROJECT_NAME}Targets"
  FILE "${PROJECT_NAME}Targets.cmake"
  NAMESPACE dhb::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}
)

add_library(${PROJECT_NAME}::${PROJECT_NAME} STATIC IMPORTED)

set_target_properties(${PROJECT_NAME}::${PROJECT_NAME} PROPERTIES
  INTERFACE_INCLUDE_DIRECTORIES "${_IMPORT_PREFIX}/include"
)

set_target_properties(dhb PROPERTIES CXX_STANDARD 17)
set_target_properties(dhb PROPERTIES POSITION_INDEPENDENT_CODE ON)
set_target_properties(dhb PROPERTIES DEBUG_POSTFIX "d")
set_target_properties(dhb PROPERTIES PUBLIC_HEADER "${public_headers}")

# ===================================
# Make a Configuration Package
# ===================================
include(CMakePackageConfigHelpers)

export(EXPORT ${PROJECT_NAME}Targets
       FILE "${CMAKE_CURRENT_BINARY_DIR}/cmake/${PROJECT_NAME}Targets.cmake"
       NAMESPACE ${PROJECT_NAME}::
)

configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/Config.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}
)

write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
  VERSION "${version}"
  COMPATIBILITY AnyNewerVersion
)

install(FILES
  "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
  "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}
)

set(version 1.0.0)

set_property(TARGET dhb PROPERTY VERSION ${version})
set_property(TARGET dhb PROPERTY SOVERSION 1)
set_property(TARGET dhb PROPERTY
  INTERFACE_dhb_MAJOR_VERSION 1)
set_property(TARGET dhb APPEND PROPERTY
  COMPATIBLE_INTERFACE_STRING dhb_MAJOR_VERSION
)

# ========================
# Test Target
# ========================
option(DHB_BUILD_TESTS "Build test target." OFF)
if (DHB_BUILD_TESTS)
  if(EXISTS "${PROJECT_SOURCE_DIR}/test/lib/Catch2")
    add_subdirectory(${PROJECT_SOURCE_DIR}/test/lib/Catch2)
  else()
    message(FATAL_ERROR
        "Missing Catch2 library in test/lib/Catch2 "
        "Please run `git submodule update --init` to fetch the submodule.")
  endif()

  if(EXISTS "${PROJECT_SOURCE_DIR}/test/graphs")
    add_definitions(-DDHB_TEST_GRAPH_DIR="${PROJECT_SOURCE_DIR}/test/graphs")
  else()
    message(FATAL_ERROR
      "Could not find test graph directory.")
  endif()

  add_executable(dhb_test
    test/graph_io.h
    test/graph_io.cpp
    test/matrix_test.cpp
    test/vec.cpp
    test/block.cpp
    test/submatrix.cpp
  )

  target_link_libraries(dhb_test PRIVATE dhb Catch2::Catch2WithMain OpenMP::OpenMP_CXX)

  install(
    TARGETS dhb_test
    RUNTIME DESTINATION bin
  )

  # Test setup. Use CTest + Catch2's CTest integration.
  include(CTest)
  include(${PROJECT_SOURCE_DIR}/test/lib/Catch2/extras/Catch.cmake)
  catch_discover_tests(dhb_test)
endif()
