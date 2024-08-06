add_library(topotoolbox
  topotoolbox.c
  fillsinks.c
  identifyflats.c
  gwdt.c
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
