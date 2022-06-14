cmake_minimum_required(VERSION 3.16)

project(PunchTable)

set(CMAKE_CXX_STANDARD 17)

if(NOT CMAKE_BUILD_TYPE  STREQUAL "Debug")
    set(CMAKE_BUILD_TYPE "Release")
    MESSAGE(STATUS "Build type Release")
else()
    MESSAGE(STATUS "Build type Debug")
endif()

set(PUNCHTABLE_BINARY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/bin)

add_subdirectory(src/vendor/catima)
add_subdirectory(src)