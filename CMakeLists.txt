cmake_minimum_required(VERSION 3.0)

set(CMAKE_CXX_STANDARD 17)

project(casadi_optimization)

find_package(PkgConfig REQUIRED)
pkg_check_modules(MATIO REQUIRED matio)

link_directories(/usr/lib/x86_64-linux-gnu/hdf5/serial)
add_executable(optimizer 3_rosenbrock.cpp)

target_include_directories(optimizer PRIVATE ${MATIO_INCLUDE_DIRS})
target_link_libraries(optimizer PRIVATE ${MATIO_LIBRARIES} casadi hdf5)

# target_include_directories(optimizer PRIVATE /usr/include/eigen3)
# target_link_directories(optimizer PRIVATE ${CMAKE_SOURCE_DIR}/casadi_api_a_test)
# target_link_libraries(optimizer PRIVATE pgd_fun)

