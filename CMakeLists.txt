cmake_minimum_required(VERSION 3.5)

project(
  topotoolbox
  VERSION 3.0.0
  LANGUAGES C CXX
)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)

include(CTest)

set(CMAKE_COMPILE_WARNING_AS_ERROR ON)
if (MSVC)
  add_compile_options(/W4 /wd4100)
else()
  add_compile_options(-Wall -Wextra -Wpedantic -Wno-unused-parameter)
endif()

add_subdirectory(src)

if (BUILD_TESTING)
  add_subdirectory(test)
endif()
