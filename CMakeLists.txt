cmake_minimum_required(VERSION 3.5)

project(
  topotoolbox
  VERSION 3.0.0
  LANGUAGES C CXX
)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)

OPTION(TT_BUILD_TESTS "Build libtopotoolbox tests" OFF)
OPTION(TT_BUILD_DOCS "Build libtopotoolbox documentation" OFF)
OPTION(TT_INSTALL "Install libtopotoolbox" OFF)

include(CTest)
include(GNUInstallDirs)

set(CMAKE_COMPILE_WARNING_AS_ERROR ON)
if (MSVC)
  add_compile_options(/W4 /wd4100)
else()
  add_compile_options(-Wall -Wextra -Wpedantic -Wno-unused-parameter -g)
endif()

add_subdirectory(src)

if (TT_BUILD_TESTS)
  add_subdirectory(test)
endif()

if (TT_BUILD_DOCS)
  add_subdirectory(docs)
endif()

if (TT_INSTALL)
  # Only set up the install targets if TT_INSTALL is set
  
  # CMake package configuration
  include(CMakePackageConfigHelpers)
  configure_package_config_file(
    cmake/TopoToolboxConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/TopoToolboxConfig.cmake
    PATH_VARS CMAKE_INSTALL_INCLUDEDIR
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}
  )

  write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/TopoToolboxConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
  )

  # Install library
  install(TARGETS topotoolbox EXPORT topotoolbox-targets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}/${PROJECT_NAME}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}/${PROJECT_NAME}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}/${PROJECT_NAME}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME})

  install(EXPORT topotoolbox-targets
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}
  )

  install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/TopoToolboxConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/TopoToolboxConfigVersion.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}
  )
endif()
