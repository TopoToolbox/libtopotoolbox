add_executable(versioninfo versioninfo.cpp utils.c)
target_link_libraries(versioninfo PRIVATE topotoolbox)
add_test(NAME versioninfo COMMAND versioninfo)
set_tests_properties(versioninfo PROPERTIES ENVIRONMENT_MODIFICATION
"PATH=path_list_prepend:$<$<BOOL:${WIN32}>:$<TARGET_FILE_DIR:topotoolbox>>")

add_executable(random_dem random_dem.cpp utils.c)
target_link_libraries(random_dem PRIVATE topotoolbox)
add_test(NAME random_dem COMMAND random_dem)
set_tests_properties(random_dem PROPERTIES ENVIRONMENT_MODIFICATION
  "PATH=path_list_prepend:$<$<BOOL:${WIN32}>:$<TARGET_FILE_DIR:topotoolbox>>")

add_executable(excesstopography excesstopography.cpp utils.c)
target_link_libraries(excesstopography PRIVATE topotoolbox)
add_test(NAME excesstopography COMMAND excesstopography)
set_tests_properties(excesstopography PROPERTIES ENVIRONMENT_MODIFICATION
  "PATH=path_list_prepend:$<$<BOOL:${WIN32}>:$<TARGET_FILE_DIR:topotoolbox>>")

add_executable(gwdt gwdt.cpp utils.c)
target_link_libraries(gwdt PRIVATE topotoolbox)
add_test(NAME gwdt COMMAND gwdt)
set_tests_properties(gwdt PROPERTIES ENVIRONMENT_MODIFICATION
  "PATH=path_list_prepend:$<$<BOOL:${WIN32}>:$<TARGET_FILE_DIR:topotoolbox>>")
