cmake_minimum_required(VERSION 3.5)
project(integral_eq)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

set(CMAKE_CXX_STANDARD 14)


add_executable(integral_eq main.cpp integral_helper.h proportion_from_k.h)
target_link_libraries(integral_eq)

