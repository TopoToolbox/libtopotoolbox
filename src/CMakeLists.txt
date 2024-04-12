add_library(topotoolbox
  topotoolbox.c
  fillsinks.c
  morphology/reconstruct.c
)

target_include_directories(
  topotoolbox
  PUBLIC
  $<BUILD_INTERFACE:${topotoolbox_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}>
)

if (BUILD_SHARED_LIBS)
else()
  target_compile_options(topotoolbox PUBLIC -DTOPOTOOLBOX_STATICLIB)
endif()

set_target_properties(topotoolbox PROPERTIES PUBLIC_HEADER ${topotoolbox_SOURCE_DIR}/include/topotoolbox.h)
