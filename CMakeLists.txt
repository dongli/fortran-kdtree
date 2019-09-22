cmake_minimum_required(VERSION 2.8)

project(kdtree LANGUAGES Fortran)

if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-none")
endif ()

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -Wall -Wextra -Warray-temporaries -Wconversion -fimplicit-none -fbacktrace -fcheck=all -ffpe-trap=zero,overflow,underflow -finit-real=nan")
else ()
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffpe-summary=none -Ofast")
endif ()

set (CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}")
include_directories(${CMAKE_BINARY_DIR})

add_subdirectory(lib/unit-test)

set(sources
  src/utils_mod.F90
  src/node_mod.F90
  src/node_placing_mod.F90
  src/qsort_mod.F90
  src/kdtree_mod.F90
)

add_library(kdtree ${sources})

add_executable(kdtree_test.exe src/kdtree_test.F90)
target_link_libraries(kdtree_test.exe kdtree fortran_unit_test)
