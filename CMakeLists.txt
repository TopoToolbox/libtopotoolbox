cmake_minimum_required(VERSION 3.5)

project(
  topotoolbox
  VERSION 3.0.0
  LANGUAGES C CXX
)

include(CTest)

add_subdirectory(src)

if (BUILD_TESTING)
  add_subdirectory(test)
endif()
