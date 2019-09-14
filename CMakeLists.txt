cmake_minimum_required(VERSION 2.8)

project(kdtree LANGUAGES Fortran)

if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffpe-summary=none -fdefault-real-8 -ffree-line-length-none -Ofast -fopenmp")
endif ()

set(sources
  src/utils_mod.F90
  src/node_mod.F90
  src/qsort_mod.F90
  src/kdtree_mod.F90
)

add_library(kdtree ${sources})

add_executable(kdtree_test.exe src/kdtree_test.F90)
target_link_libraries(kdtree_test.exe kdtree)
