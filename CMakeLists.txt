cmake_minimum_required(VERSION 3.6)
project(bh)

set(CMAKE_CXX_STANDARD 11)

find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

set(CMAKE_CXX_COMPILER mpicxx)

set(SOURCE_FILES main.cpp body.cpp body.h quad.cpp quad.h tree_buffer.cpp tree_buffer.h)

add_executable(bh ${SOURCE_FILES})