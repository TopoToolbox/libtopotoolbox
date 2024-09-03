extern "C" {
#include "topotoolbox.h"
}

#include <iostream>

int main(int argc, char *argv[]) {
  if (has_topotoolbox()) {
    std::cout << "topotoolbox v" << TOPOTOOLBOX_VERSION_MAJOR << "."
              << TOPOTOOLBOX_VERSION_MINOR << "." << TOPOTOOLBOX_VERSION_PATCH
              << std::endl;

    std::cout << "openmp v" << TOPOTOOLBOX_OPENMP_VERSION << std::endl;
  }
}
