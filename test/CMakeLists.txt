add_executable(versioninfo versioninfo.cpp)
target_link_libraries(versioninfo PRIVATE topotoolbox)
add_test(NAME versioninfo COMMAND versioninfo)
set_tests_properties(versioninfo PROPERTIES ENVIRONMENT_MODIFICATION
"PATH=path_list_prepend:$<$<BOOL:${WIN32}>:$<TARGET_FILE_DIR:topotoolbox>>")

add_executable(fillsinks fillsinks.cpp)
target_link_libraries(fillsinks PRIVATE topotoolbox)
add_test(NAME fillsinks COMMAND fillsinks)
set_tests_properties(fillsinks PROPERTIES ENVIRONMENT_MODIFICATION
"PATH=path_list_prepend:$<$<BOOL:${WIN32}>:$<TARGET_FILE_DIR:topotoolbox>>")
