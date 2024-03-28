add_library(topotoolbox SHARED
  topotoolbox.c
  fillsinks.c
  morphology/reconstruct.c
)

target_include_directories(topotoolbox PUBLIC ../include)
