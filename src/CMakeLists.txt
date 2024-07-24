add_library(topotoolbox
  topotoolbox.c
  fillsinks.c
  identifyflats.c
  gwdt.c
  morphology/reconstruct.c
  priority_queue.c
  excesstopography.c
  graphflood/graphflood.c
)

target_include_directories(
  topotoolbox
  PUBLIC
  $<BUILD_INTERFACE:${topotoolbox_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}>
)

find_library(MATH_LIBRARY m)
if(MATH_LIBRARY)
    target_link_libraries(topotoolbox PUBLIC ${MATH_LIBRARY})
endif()

if (BUILD_SHARED_LIBS)
else()
  target_compile_options(topotoolbox PUBLIC -DTOPOTOOLBOX_STATICLIB)
endif()

set_target_properties(topotoolbox PROPERTIES PUBLIC_HEADER ${topotoolbox_SOURCE_DIR}/include/topotoolbox.h)
