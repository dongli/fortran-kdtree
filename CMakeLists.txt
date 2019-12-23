cmake_minimum_required(VERSION 2.8)

project(kdtree LANGUAGES Fortran)

if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-none")
  if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -Wall -Wextra -Warray-temporaries -Wconversion -fimplicit-none -fbacktrace -fcheck=all -ffpe-trap=zero,overflow,underflow -finit-real=nan")
  else ()
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffpe-summary=none -Ofast")
  endif ()
endif ()

set (CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}")
include_directories(${CMAKE_BINARY_DIR})

if (EXISTS ${PROJECT_SOURCE_DIR}/lib/unit-test/CMakeLists.txt)
  add_subdirectory(lib/unit-test)
  set(RUN_TEST ON)
endif ()

set(sources
  src/utils_mod.F90
  src/node_mod.F90
  src/node_placing_mod.F90
  src/qsort_mod.F90
  src/kdtree_mod.F90
  src/kdtree.F90
)

add_library(kdtree ${sources})

if (RUN_TEST)
  add_executable(kdtree_test.exe src/kdtree_test.F90)
  target_link_libraries(kdtree_test.exe kdtree fortran_unit_test)
endif ()
