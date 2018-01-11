cmake_minimum_required(VERSION 3.6)
project(btl)

set(CMAKE_CXX_STANDARD 11)

find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

set(CMAKE_CXX_COMPILER mpicxx)

set(SOURCE_FILES main.cpp)

add_executable(btl ${SOURCE_FILES})