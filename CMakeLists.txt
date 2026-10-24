# morphe: mathematical structures and their equivalences.
#
# Every module is a target morphe::<module>. morphe depends on no other library of the
# project, not even for its checks: they print the shared report through morphe/check.hpp
# (decision 0004). Built on its own, morphe builds its checks (MORPHE_CHECKS); added to a
# larger project, it builds them only when asked.

cmake_minimum_required(VERSION 3.20)
project(morphe VERSION 0.1.0 LANGUAGES CXX)

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  set(MORPHE_TOP_LEVEL ON)
else()
  set(MORPHE_TOP_LEVEL OFF)
endif()
option(MORPHE_CHECKS "Build morphe's checks" ${MORPHE_TOP_LEVEL})

if(NOT DEFINED CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 20)
  set(CMAKE_CXX_STANDARD_REQUIRED ON)
  set(CMAKE_CXX_EXTENSIONS OFF)
endif()
if(MORPHE_CHECKS)
  enable_testing()
endif()

# The report of the checks, shared by every module's check program.
add_library(morphe_check INTERFACE)
add_library(morphe::check ALIAS morphe_check)
target_include_directories(morphe_check INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_features(morphe_check INTERFACE cxx_std_20)

# morphe_add_check(<name> SOURCES <files> LIBRARIES <targets>): a check program, run by CTest.
function(morphe_add_check name)
  cmake_parse_arguments(ARG "" "" "SOURCES;LIBRARIES" ${ARGN})
  add_executable(${name} ${ARG_SOURCES})
  target_link_libraries(${name} PRIVATE morphe::check ${ARG_LIBRARIES})
  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(${name} PRIVATE -Wall -Wextra)
  endif()
  add_test(NAME ${name} COMMAND ${name})
endfunction()

# The modules in their order (docs/architecture.md §2).
add_subdirectory(arithmetic)
add_subdirectory(algebra)
add_subdirectory(combinatorics)
add_subdirectory(topology)
add_subdirectory(formal)
