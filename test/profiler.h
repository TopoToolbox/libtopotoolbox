#ifndef TT_PROFILER_H
#define TT_PROFILER_H

/*
  A basic profiler for libtopotoolbox.

  Usage:
  1. `#include "profiler.h"` in your test executable.
  2. Create a profiler with `Profiler prof;`. It is easiest to create a
     global variable so that everything can access it, but you can also pass it
     by reference into the functions that need it.
  3. Put `ProfileBlock(prof,label);` at the start of each block scope
     in the test executable file you want to time. `label` should be a string
     literal that labels the block. You can use `ProfileFunction(prof)` to
     profile a function, in which case the label is the name of the function.
  4. Call `prof.report()` when you are finished to print the profiling results.

  The output is a JSON object with a single field "blocks" that
  points to an array of objects of the form

  `{"label" : label of block or name of function,
    "calls": number of calls of that function during the test,
    "time": average time per call in milliseconds
    }`

  Known Limitations:

  - Recursive and nested function calls can lead to misleading results.
 */

#include <cstdint>
#include <ctime>
#include <iostream>
#include <string>
#include <unordered_map>

#define ProfileBlock(Profiler, Label) \
  ProfileZone block##__LINE__(Profiler, Label);
#define ProfileFunction(Profiler) \
  ProfileZone block##__LINE__(Profiler, __func__);

struct ProfileStats {
  uint64_t elapsed;
  uint64_t count;
};

struct Profiler {
  void report() {
    std::cout << "{\"blocks\": [" << std::endl;
    int count = 0;
    for (auto iter = anchors.begin(); iter != anchors.end(); ++iter) {
      if (iter->second.count > 0) {
        if (count != 0) {
          std::cout << "," << std::endl;
        }
        double anchor_ms = 1000.0 * iter->second.elapsed /
                           (double)CLOCKS_PER_SEC / iter->second.count;

        std::cout << "{\"label\": \"" << iter->first << "\"," << std::endl;
        std::cout << "\"calls\": " << iter->second.count << "," << std::endl;
        std::cout << "\"time\": " << anchor_ms << "}" << std::endl;
        count++;
      }
    }
    std::cout << "]}" << std::endl;
  }

  ProfileStats &operator[](std::string label) { return anchors[label]; }

 private:
  std::unordered_map<std::string, ProfileStats> anchors;
};

struct ProfileZone {
  ProfileZone(Profiler &prof_, const char *label_)
      : label(label_), prof(prof_) {
    start = std::clock();
  }

  ~ProfileZone(void) {
    uint64_t elapsed = std::clock() - start;

    ProfileStats &anchor = prof[label];
    anchor.elapsed += elapsed;
    anchor.count++;
  }

 private:
  const char *label;
  uint64_t start;
  Profiler &prof;
};

#endif  // TT_PROFILER_H