cmake_minimum_required(VERSION 3.22)
project(6DOF_simplified)

file(GLOB_RECURSE SOURCE_FILES src/*.cpp)
add_executable(main parafoil_simulator.cpp ${SOURCE_FILES})

target_compile_features(main PRIVATE cxx_std_17)
target_include_directories(main PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
#boost
find_package(Boost REQUIRED)
target_link_libraries(main ${Boost_LIBRARIES})
include_directories(${Boost_INCLUDE_DIR})

