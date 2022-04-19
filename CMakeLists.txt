cmake_minimum_required(VERSION 3.0)
project(grizz_ins C)

set(CMAKE_C_STANDARD 99)

link_libraries(m)

add_executable(grizz_ins main.c inertial.h inertial.c inertial_math.c inertial_math.h platform.h)
