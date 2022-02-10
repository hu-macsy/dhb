# Set minimum version of CMake.
cmake_minimum_required(VERSION 3.11.4)

# Set project name and version
project(dhb VERSION 1.0.0)

include(GNUInstallDirs)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Configure the option whether you want to use the systems own allocater
# WITH_DHB_SYSTEM_ALLOCATOR=ON
# or our own written allocator strategy
# WITH_DHB_SYSTEM_ALLOCATOR=OFF
option(WITH_DHB_SYSTEM_ALLOCATOR "Use DHB System Allocator" OFF)
option(WITH_DHB_64BIT_IDS "Use 64 bit IDs." OFF)

# different ways of sorting edges when inserting in concurrently
set(WITH_DHB_SCATTER "twophase" CACHE STRING "Parallelization method")

add_library(dhb STATIC)

set(public_headers
  include/dhb/batcher.h
  include/dhb/block.h
  include/dhb/buckets.h
  include/dhb/dynamic_hashed_blocks.h
  include/dhb/graph.h
  include/dhb/integer_log2.h
  include/dhb/vec.h
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
if (WITH_DHB_SYSTEM_ALLOCATOR)
  MESSAGE(STATUS "Using DHB System Allocator.")
  target_compile_definitions(dhb PUBLIC -DDHB_SYSTEM_ALLOCATOR)
endif()

if (WITH_DHB_64BIT_IDS)
  MESSAGE(STATUS "Using 64 bit IDs.")
  target_compile_definitions(dhb PUBLIC -DDHB_64BIT_IDS)
else()
  MESSAGE(STATUS "Using 32 bit IDs.")
endif()

if (WITH_DHB_SCATTER STREQUAL "twophase")
	MESSAGE(STATUS "Using scatter in two phases.")
	target_compile_definitions(dhb PUBLIC -DDHB_SCATTER_TWOPHASE)
elseif (WITH_DHB_SCATTER STREQUAL "counting")
	MESSAGE(STATUS "Using scatter counting.")
	target_compile_definitions(dhb PUBLIC -DDHB_SCATTER_COUNTING)
elseif (WITH_DHB_SCATTER STREQUAL "sorting")
	MESSAGE(STATUS "Using scatter sorting.")
	target_compile_definitions(dhb PUBLIC -DDHB_SCATTER_SORTING)
elseif (WITH_DHB_SCATTER STREQUAL "darts")
	MESSAGE(STATUS "Using scatter darts.")
	target_compile_definitions(dhb PUBLIC -DDHB_SCATTER_DARTS)
else()
	message(FATAL_ERROR "Illegal setting for WITH_DHB_SCATTER")
endif()

set_property(TARGET dhb PROPERTY POSITION_INDEPENDENT_CODE ON)

find_package(OpenMP)

target_link_libraries(${PROJECT_NAME} PRIVATE OpenMP::OpenMP_CXX)

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

set_target_properties(${PROJECT_NAME} PROPERTIES PUBLIC_HEADER "${public_headers}")
set_target_properties(${PROJECT_NAME} PROPERTIES DEBUG_POSTFIX "d")

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
option(DHB_TEST "Build test target." OFF)
if (DHB_TEST)
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
    test/matrix.cpp
    test/vec.cpp
    test/block.cpp
  )

  target_link_libraries(dhb_test PRIVATE dhb Catch2::Catch2WithMain OpenMP::OpenMP_CXX stdc++fs)

  install(
    TARGETS dhb_test
    RUNTIME DESTINATION bin
  )

  # Test setup. Use CTest + Catch2's CTest integration.
  include(CTest)
  include(${PROJECT_SOURCE_DIR}/test/lib/Catch2/extras/Catch.cmake)
  catch_discover_tests(dhb_test)
endif()
