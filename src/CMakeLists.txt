# Library declaration
#
# Every file that needs to be compiled should be listed here.
add_library(topotoolbox
  topotoolbox.c
  fillsinks.c
  identifyflats.c
  gwdt.c
  gradient8.c
  morphology/reconstruct.c
  priority_queue.c
  excesstopography.c
  graphflood/gf_utils.c
  graphflood/sfgraph.c
  graphflood/priority_flood_standalone.c
  graphflood/gf_flowacc.c
  graphflood/graphflood.c
  flow_routing.c
  flow_accumulation.c
)

# Define the include directory
#
# When we are building the library, look in include/
#
# When we install the library, look in the appropriate system include
# directory
target_include_directories(
  topotoolbox
  PUBLIC
  $<BUILD_INTERFACE:${topotoolbox_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}>
)

# Link libm if necessary
#
# Not every platform has a separate libm that needs to be linked, this
# conditional ensures that we only try to link libm if it exists separately

find_library(MATH_LIBRARY m)
if(MATH_LIBRARY)
    target_link_libraries(topotoolbox PUBLIC ${MATH_LIBRARY})
endif()

# Link OpenMP if found
#
# This conditional ensures that we only try to link OpenMP if it is found
# on the system. The OpenMP_C_FOUND variable is set by the find_package
# command and indicates whether the OpenMP C library is available.

find_package(OpenMP)
if (OpenMP_C_FOUND)
  target_compile_options(topotoolbox PUBLIC -DTOPOTOOLBOX_OPENMP_VERSION=${OpenMP_C_VERSION_MAJOR}${OpenMP_C_VERSION_MINOR})
  target_link_libraries(topotoolbox PUBLIC OpenMP::OpenMP_C)
endif()

# Set TOPOTOOLBOX_STATICLIB if we are *not* building a shared library
#
# TOPOTOOLBOX_STATICLIB is a macro defined in include/topotoolbox.h
# that determines the appropriate TOPOTOOLBOX_API qualifier when
# libtopotoolbox is built as a static library.

if (BUILD_SHARED_LIBS)
else()
  target_compile_options(topotoolbox PUBLIC -DTOPOTOOLBOX_STATICLIB)
endif()

# Define include/topotoolbox.h as a PUBLIC_HEADER
#
# This is used in the root CMakeLists.txt to install topotoolbox.h to
# the appropriate location.
set_target_properties(topotoolbox PROPERTIES PUBLIC_HEADER ${topotoolbox_SOURCE_DIR}/include/topotoolbox.h)
