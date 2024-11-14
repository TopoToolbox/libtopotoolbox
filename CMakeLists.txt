cmake_minimum_required(VERSION 3.12)

project(
  topotoolbox
  VERSION 3.0.0
  LANGUAGES C CXX
)

# libtopotoolbox/cmake contains additional CMake modules necessary for
# our build process.
#
# If you need build system functionality that is not built into CMake, add
# the relevant CMake module to this directory.
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

# Language standards and compiler settings
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)

# TopoToolbox specific build options
#
# These should start with the prefix `TT_`
OPTION(TT_BUILD_TESTS "Build libtopotoolbox tests" OFF)
OPTION(TT_BUILD_DOCS "Build libtopotoolbox documentation" OFF)
OPTION(TT_INSTALL "Install libtopotoolbox" OFF)


# Compiler warnings
set(CMAKE_COMPILE_WARNING_AS_ERROR ON)
if (MSVC)
  add_compile_options(/W4 /wd4100)
else()
  add_compile_options(-Wall -Wextra -Wpedantic -Wno-unused-parameter)
endif()

# Library source code is in the src/ directory
# See src/CMakeLists.txt for the library configuration
add_subdirectory(src)

# Tests are in the test/ directory
# See test/CMakeLists.txt for the test configuration
include(CTest)

if (TT_BUILD_TESTS)
  add_subdirectory(test)
endif()

# Documentation is built from the docs/ directory
# See docs/CMakeLists.txt for the documentation configuration
if (TT_BUILD_DOCS)
  add_subdirectory(docs)
endif()

# Set up the install targets
include(GNUInstallDirs)

if (TT_INSTALL)
  # Only set up the install targets if TT_INSTALL is set
  
  # CMake package configuration
  #
  # Generate TopoToolboxConfig.cmake, which find_package()
  # uses to locate an installed version of libtopotoolbox.
  #
  # This file will be installed in something like
  # $PREFIX/lib/cmake/topotoolbox, depending on the OS.
  include(CMakePackageConfigHelpers)
  configure_package_config_file(
    cmake/TopoToolboxConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/TopoToolboxConfig.cmake
    PATH_VARS CMAKE_INSTALL_INCLUDEDIR
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}
  )

  # Generate the TopoToolboxConfigVersion.cmake, which find_package()
  # uses to ensure that the requested version is installed.
  write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/TopoToolboxConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
  )

  # Install the compiled library and header file
  install(TARGETS topotoolbox EXPORT topotoolbox-targets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}/${PROJECT_NAME}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}/${PROJECT_NAME}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}/${PROJECT_NAME}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME})

  install(EXPORT topotoolbox-targets
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}
  )

  # Install the CMake package configuration files in the appropriate
  # location
  install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/TopoToolboxConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/TopoToolboxConfigVersion.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}
  )
endif()
