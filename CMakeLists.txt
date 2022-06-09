cmake_minimum_required(VERSION 3.16)

project(PunchTable)

set(CXX_STANDARD 17)

set(PUNCHTABLE_BINARY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/bin)

add_subdirectory(src/vendor/catima)
add_subdirectory(src)