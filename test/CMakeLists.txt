# Set TT_SANITIZE=ON to build the tests with AddressSanitizer.
#
# We run the tests with AddressSanitizer on CI to catch certain memory
# access errors, but AddressSanitizer adds significant overhead to
# programs, so it is not enabled by default.
OPTION(TT_SANITIZE "Build with AddressSanitizer" OFF)

# Test declarations should generally follow the pattern below, which
# is described in more detail for the versioninfo test.

# TEST : versioninfo
#
# This test prints out the defined version numbers. It is mainly used
# to ensure that everything is linked properly.

# add_executable defines the executable program that CTest will run.
# Include any necessary files that must be compiled in this command.
add_executable(versioninfo versioninfo.cpp)

# Link the executable to topotoolbox.
target_link_libraries(versioninfo PRIVATE topotoolbox)

# Add the executable to CTest so that it will run it.
add_test(NAME versioninfo COMMAND versioninfo)

# This property is required on Windows to ensure that the shared
# library is found when running the test.
set_tests_properties(versioninfo PROPERTIES ENVIRONMENT_MODIFICATION
                     "PATH=path_list_prepend:$<$<BOOL:${WIN32}>:$<TARGET_FILE_DIR:topotoolbox>>")

# TEST : random_dem
#
# Runs libtopotoolbox functions on randomly generated DEMs.
add_executable(random_dem random_dem.cpp utils.c utils.h profiler.h)
if(TT_SANITIZE AND NOT MSVC)
  # Set up the AddressSanitizer flags.
  target_compile_options(random_dem PRIVATE "$<$<CONFIG:DEBUG>:-fsanitize=address>")
  target_link_options(random_dem PRIVATE "$<$<CONFIG:DEBUG>:-fsanitize=address>")
endif()
target_link_libraries(random_dem PRIVATE topotoolbox)
add_test(NAME random_dem COMMAND random_dem)
set_tests_properties(random_dem PROPERTIES ENVIRONMENT_MODIFICATION
                     "PATH=path_list_prepend:$<$<BOOL:${WIN32}>:$<TARGET_FILE_DIR:topotoolbox>>")

# TEST : excesstopography
#
# Runs the excesstopography functions on randomly generated DEMs.
add_executable(excesstopography excesstopography.cpp utils.c utils.h)
if(TT_SANITIZE AND NOT MSVC)
  # Set up the AddressSanitizer flags.
  target_compile_options(excesstopography PRIVATE "$<$<CONFIG:DEBUG>:-fsanitize=address>")
  target_link_options(excesstopography PRIVATE "$<$<CONFIG:DEBUG>:-fsanitize=address>")
endif()
target_link_libraries(excesstopography PRIVATE topotoolbox)
add_test(NAME excesstopography COMMAND excesstopography)
set_tests_properties(excesstopography PROPERTIES ENVIRONMENT_MODIFICATION
                     "PATH=path_list_prepend:$<$<BOOL:${WIN32}>:$<TARGET_FILE_DIR:topotoolbox>>")

# TEST : filters
#
# This test runs various property tests on value filters and morphological
# operations on randomly generated DEMs.
add_executable(filters filters.cpp utils.c)
target_link_libraries(filters PRIVATE topotoolbox)
add_test(NAME filters COMMAND filters)
set_tests_properties(filters PROPERTIES ENVIRONMENT_MODIFICATION
                     "PATH=path_list_prepend:$<$<BOOL:${WIN32}>:$<TARGET_FILE_DIR:topotoolbox>>")

# TEST : swaths
#
# Tests swath_frontier_distance_map, swath_boundary_dijkstra,
# voronoi_ridge_to_centreline, swath_longitudinal, and
# swath_get_point_pixels on a straight horizontal track with
# analytically known distances.
add_executable(swaths swaths.cpp)
if(TT_SANITIZE AND NOT MSVC)
  target_compile_options(swaths PRIVATE "$<$<CONFIG:DEBUG>:-fsanitize=address>")
  target_link_options(swaths PRIVATE "$<$<CONFIG:DEBUG>:-fsanitize=address>")
endif()
target_link_libraries(swaths PRIVATE topotoolbox)
add_test(NAME swaths COMMAND swaths)
set_tests_properties(swaths PROPERTIES ENVIRONMENT_MODIFICATION
                     "PATH=path_list_prepend:$<$<BOOL:${WIN32}>:$<TARGET_FILE_DIR:topotoolbox>>")

# TEST : polyline
#
# Tests rasterize_path, thin_rasterised_line_to_D8, and simplify_line.
add_executable(polyline polyline.cpp)
if(TT_SANITIZE AND NOT MSVC)
  target_compile_options(polyline PRIVATE "$<$<CONFIG:DEBUG>:-fsanitize=address>")
  target_link_options(polyline PRIVATE "$<$<CONFIG:DEBUG>:-fsanitize=address>")
endif()
target_link_libraries(polyline PRIVATE topotoolbox)
add_test(NAME polyline COMMAND polyline)
set_tests_properties(polyline PROPERTIES ENVIRONMENT_MODIFICATION
                     "PATH=path_list_prepend:$<$<BOOL:${WIN32}>:$<TARGET_FILE_DIR:topotoolbox>>")

# TEST : snapshots
#
# Runs the snapshot tests from available snapshot data.
OPTION(TT_DOWNLOAD_SNAPSHOTS
       "Download snapshot test data if it is not already present in test/snapshots/data" OFF)

set(TT_SNAPSHOT_VERSION v1.7.0)

# This is the SHA256 hash of the above version's snapshot_data.tar.gz archive.
set(TT_SNAPSHOT_HASH 77024446fa7764ce68e89d0e8cc32adf73df1953d4100fa4c861321fbc6c0c47)

set(TT_SNAPSHOT_DATA_PATH ${CMAKE_CURRENT_SOURCE_DIR}/snapshots/data)

include(FetchContent)

if(EXISTS ${TT_SNAPSHOT_DATA_PATH})
  message(STATUS "Found existing snapshot data, copying it into the build directory")
  file(COPY ${TT_SNAPSHOT_DATA_PATH} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/snapshots/)
  set(TT_SNAPSHOT_DIR ${CMAKE_CURRENT_BINARY_DIR}/snapshots/data)
  message(STATUS "Snapshot data directory: ${TT_SNAPSHOT_DIR}")
elseif(TT_DOWNLOAD_SNAPSHOTS)
  # If the snapshot data does not yet exist at test/snapshot/data, but
  # we still want to download the snapshots, we do it with FetchContent.
  message(STATUS "Downloading snapshot data")
  FetchContent_Declare(
    snapshot_data
    URL https://github.com/TopoToolbox/snapshot_data/releases/download/${TT_SNAPSHOT_VERSION}/snapshot_data.tar.gz
    URL_HASH SHA256=${TT_SNAPSHOT_HASH}
    SOURCE_DIR ${TT_SNAPSHOT_DATA_PATH}
    DOWNLOAD_EXTRACT_TIMESTAMP false)
  FetchContent_MakeAvailable(snapshot_data)
  set(TT_SNAPSHOT_DIR ${snapshot_data_SOURCE_DIR})
  message(STATUS "Snapshot data directory: ${TT_SNAPSHOT_DIR}")
else()
  # If we did not download the snapshot data, and the snapshot data
  # was not provided, we unset TT_SNAPSHOT_DIR so that we do not build
  # the snapshot test.
  unset(TT_SNAPSHOT_DIR)
  message(STATUS "Could not find snapshot data. Snapshot tests will not be run")
endif()

include(FindGDAL)
if(GDAL_FOUND AND TT_SNAPSHOT_DIR)
  add_executable(snapshot snapshot.cpp utils.c utils.h profiler.h)
  if(TT_SANITIZE AND NOT MSVC)
    # Set up the AddressSanitizer flags.
    target_compile_options(snapshot PRIVATE "$<$<CONFIG:DEBUG>:-fsanitize=address>")
    target_link_options(snapshot PRIVATE "$<$<CONFIG:DEBUG>:-fsanitize=address>")
  endif()
  target_link_libraries(snapshot PRIVATE topotoolbox GDAL::GDAL)
  add_test(NAME snapshot COMMAND snapshot ${TT_SNAPSHOT_DIR})
  set_tests_properties(snapshot PROPERTIES ENVIRONMENT_MODIFICATION
                       "PATH=path_list_prepend:$<$<BOOL:${WIN32}>:$<TARGET_FILE_DIR:topotoolbox>>")
endif()

# TEST: clang-format
#
# Runs clang format to verify formatting of code

OPTION(TT_FORMAT_CHECK "Check code formatting with clang-format" ON)

find_program(CLANG_FORMAT_EXE clang-format)
if (TT_FORMAT_CHECK AND CLANG_FORMAT_EXE)
  message(STATUS "Found clang-format at ${CLANG_FORMAT_EXE}")

  set(FORMAT_TARGETS
    topotoolbox
    versioninfo
    random_dem
    excesstopography
    filters
    swaths
    polyline)

  if (TARGET snapshot)
    list(APPEND FORMAT_TARGETS snapshot)
  endif()

  set(ALL_SOURCES)

  foreach(FORMAT_TARGET ${FORMAT_TARGETS})
    get_target_property(TARGET_SOURCES ${FORMAT_TARGET} SOURCES)
    get_target_property(TARGET_SOURCE_DIR ${FORMAT_TARGET} SOURCE_DIR)
    foreach (source ${TARGET_SOURCES})
      if(IS_ABSOLUTE ${source})
	list(APPEND ALL_SOURCES ${source})
      else()
	list(APPEND ALL_SOURCES ${TARGET_SOURCE_DIR}/${source})
      endif()
    endforeach()
  endforeach()
  list(REMOVE_DUPLICATES ALL_SOURCES)

  add_test(NAME clang-format
    COMMAND ${CLANG_FORMAT_EXE} -n -Werror ${ALL_SOURCES} ${CMAKE_SOURCE_DIR}/include/topotoolbox.h
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  )
endif()
