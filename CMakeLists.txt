cmake_minimum_required(VERSION 2.8)

project(Mat C CXX)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
include_directories(${CMAKE_CURRENT_LIST_DIR}/src)
file(GLOB MAT_SRC ${CMAKE_CURRENT_LIST_DIR}/src/*.h
        ${CMAKE_CURRENT_LIST_DIR}/src/*.cpp)
set(MAT_COMPILE_CODE ${MAT_SRC})
add_executable(Mat ${MAT_COMPILE_CODE})
