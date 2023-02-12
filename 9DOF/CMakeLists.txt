cmake_minimum_required(VERSION 3.22)
project(9DOF)

file(GLOB_RECURSE SOURCE_FILES src/*.cpp)
add_executable(9dof 9dof_simulator.cpp ${SOURCE_FILES})

target_compile_features(9dof PRIVATE cxx_std_17)
target_include_directories(9dof PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)