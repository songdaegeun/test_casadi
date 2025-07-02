cmake_minimum_required(VERSION 3.0)
project(casadi_optimization)

add_executable(optimizer 4_rosenbrock.cpp)
target_link_libraries(optimizer casadi)

# message(STATUS "CASADI_INCLUDE_DIRS = ${CASADI_INCLUDE_DIRS}")
# message(STATUS "CASADI_LIBRARIES = ${CASADI_LIBRARIES}")

