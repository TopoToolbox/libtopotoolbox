# Set TT_SANITIZE=ON to build the tests with AddressSanitizer
#
# We run the tests with AddressSanitizer on CI to catch certain memory
# access errors, but AddressSanitizer adds significant overhead to
# programs, so it is not enabled by default.
OPTION(TT_SANITIZE "Build with AddressSanitizer" OFF)

# Test declarations should generally follow the pattern below, which
# is described in more detail for the versioninfo test

# TEST: versioninfo
#
# This test prints out the defined version numbers. It is mainly used
# to ensure that everything is linked properly.

# add_executable defines the executable program that CTest will
# run. Include any necessary files that must be compiled in this command.
add_executable(versioninfo versioninfo.cpp utils.c)

# Link the executable to topotoolbox.
target_link_libraries(versioninfo PRIVATE topotoolbox)

# Add the executable to CTest so that it will run it.
add_test(NAME versioninfo COMMAND versioninfo)

# This property is required on Windows to ensure that the shared
# library is found when running the test.
set_tests_properties(versioninfo PROPERTIES ENVIRONMENT_MODIFICATION
"PATH=path_list_prepend:$<$<BOOL:${WIN32}>:$<TARGET_FILE_DIR:topotoolbox>>")

# TEST: random_dem
#
# Runs libtopotoolbox functions on randomly generated DEMs
add_executable(random_dem random_dem.cpp utils.c)
if (TT_SANITIZE AND NOT MSVC)
  # Set up the AddressSanitizer flags
  target_compile_options(random_dem PRIVATE "$<$<CONFIG:DEBUG>:-fsanitize=address>")
  target_link_options(random_dem PRIVATE "$<$<CONFIG:DEBUG>:-fsanitize=address>")
endif()
target_link_libraries(random_dem PRIVATE topotoolbox)
add_test(NAME random_dem COMMAND random_dem)
set_tests_properties(random_dem PROPERTIES ENVIRONMENT_MODIFICATION
  "PATH=path_list_prepend:$<$<BOOL:${WIN32}>:$<TARGET_FILE_DIR:topotoolbox>>")

# TEST: excesstopography
#
# Runs the excesstopography functions on randomly generated DEMs.
add_executable(excesstopography excesstopography.cpp utils.c)
if (TT_SANITIZE AND NOT MSVC)
  # Set up the AddressSanitizer flags
  target_compile_options(random_dem PRIVATE "$<$<CONFIG:DEBUG>:-fsanitize=address>")
  target_link_options(random_dem PRIVATE "$<$<CONFIG:DEBUG>:-fsanitize=address>")
endif()
target_link_libraries(excesstopography PRIVATE topotoolbox)
add_test(NAME excesstopography COMMAND excesstopography)
set_tests_properties(excesstopography PROPERTIES ENVIRONMENT_MODIFICATION
  "PATH=path_list_prepend:$<$<BOOL:${WIN32}>:$<TARGET_FILE_DIR:topotoolbox>>")

# TEST: snapshots
#
# Runs the snapshot tests from available snapshot data
include(FindGDAL)

if (GDAL_FOUND)
  add_executable(snapshot snapshot.cpp utils.c)
  if (TT_SANITIZE AND NOT MSVC)
    # Set up the AddressSanitizer flags
    target_compile_options(random_dem PRIVATE "$<$<CONFIG:DEBUG>:-fsanitize=address>")
    target_link_options(random_dem PRIVATE "$<$<CONFIG:DEBUG>:-fsanitize=address>")
  endif()
  target_link_libraries(snapshot PRIVATE topotoolbox GDAL::GDAL)

  # Copy test snapshots to output
  add_custom_target(copy-test-files ALL
    ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/snapshots ${CMAKE_CURRENT_BINARY_DIR}/snapshots)

  add_test(NAME snapshot COMMAND snapshot snapshots/data)
  set_tests_properties(snapshot PROPERTIES ENVIRONMENT_MODIFICATION
    "PATH=path_list_prepend:$<$<BOOL:${WIN32}>:$<TARGET_FILE_DIR:topotoolbox>>")
endif()
