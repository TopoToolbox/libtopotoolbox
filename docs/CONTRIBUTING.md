# Contribution Guidelines

Thanks for your interest in contributing to TopoToolbox!

This repository is for libtopotoolbox, the C/C++ library that powers
TopoToolbox v3. If you are interested in implementing high-performance
algorithms for terrain analysis, this is the place for you. If you are
more interested in analysis and visualization tools built on top of
libtopotoolbox, you might want to check out the bindings to
libtopotoolbox in higher-level languages such as Python
([pytopotoolbox](https://github/TopoToolbox/pytopotoolbox)) and MATLAB
([topotoolbox](https://github/TopoToolbox/topotoolbox)). If you are
interested mainly in using TopoToolbox for your work, studies or
research, you might find the documentation at
<https://topotoolbox.github.io> more useful.

We use git as a version control system and GitHub's Issues and Pull
Requests to manage contributions to libtopotoolbox. If you need some
help using git or GitHub, check out GitHub's documentation: 
[About GitHub and Git](https://docs.github.com/en/get-started/start-your-journey/about-github-and-git). 

You will need a GitHub account to contribute through Issues or Pull
Requests. If you are unable to use GitHub, patches may also be [sent
by email](https://git-send-email.io/) to
william.kearney@uni-potsdam.de.

## Opening an issue

If you would like to report a bug or request a feature you can [create
an
issue](https://github.com/TopoToolbox/libtopotoolbox/issues/new). Search
the existing issues, including closed issues, to see if your problem
has been discussed before. 

If you are reporting a bug, please provide

- A full error message
- A minimal working example of code that triggers the error
- Details such as the operating system and the C compiler used to
  compile libtopotoolbox, if available.

## Contributing via pull requests

If you would like to contribute code, tests, or documentation to
libtopotoolbox, it is a good idea to first check the issue tracker and
discuss your contribution in a new or existing issue there. It is not
necessary to get approval before making a contribution, but checking
in can help make sure that you are not actively duplicating effort or
spending time working on a feature that is unlikely to be merged.

Create your own
[fork](https://github.com/TopoToolbox/libtopotoolbox/fork) of the
repository. Make the necessary changes to a copy of the repository on
your local machine, commit them, and push them to your fork on
GitHub. Then open a pull request to the TopoToolbox/libtopotoolbox
repository. One of the TopoToolbox maintainers will review your pull
request and let you know if any changes need to be made before your
contribution can be accepted. Requested changes can be pushed to your
fork and will automatically show up in the pull request.

We typically will merge pull requests using GitHub's "Squash and
Merge" option, which combines all of the commits in a single pull
request into one large commit and then makes a merge commit to the
`main` branch of libtopotoolbox. Note that once your pull request is
merged in this way, your fork will have diverged with the upstream
TopoToolbox repository, and you will have to merge or rebase the
upstream branch into your local copy.

GitHub Actions will run a series of automated tests and formatting
checks on your pull request and let you know if any of them fail.

### Licensing and Developer Certificate of Origin

libtopotoolbox is an open source project licensed under the GNU Public
License v3.0 (see [LICENSE](../LICENSE)). We use the [Developer
Certificate of Origin](https://developercertificate.org/) to ensure
that incoming contributions are correctly attributed and licensed.

```
Developer Certificate of Origin
Version 1.1

Copyright (C) 2004, 2006 The Linux Foundation and its contributors.

Everyone is permitted to copy and distribute verbatim copies of this
license document, but changing it is not allowed.


Developer's Certificate of Origin 1.1

By making a contribution to this project, I certify that:

(a) The contribution was created in whole or in part by me and I
    have the right to submit it under the open source license
    indicated in the file; or

(b) The contribution is based upon previous work that, to the best
    of my knowledge, is covered under an appropriate open source
    license and I have the right under that license to submit that
    work with modifications, whether created in whole or in part
    by me, under the same open source license (unless I am
    permitted to submit under a different license), as indicated
    in the file; or

(c) The contribution was provided directly to me by some other
    person who certified (a), (b) or (c) and I have not modified
    it.

(d) I understand and agree that this project and the contribution
    are public and that a record of the contribution (including all
    personal information I submit with it, including my sign-off) is
    maintained indefinitely and may be redistributed consistent with
    this project or the open source license(s) involved.
```

To indicate your agreement with the DCO, please add a line in your git
commits that looks like

```
Signed-off-by: William Kearney <william.kearney@uni-potsdam.de>
```

-- substituting your name and email address. Git can also add this
statement automatically if you use the `-s` flag on the command line
of `git commit`.

## Developing libtopotoolbox

### Building and testing the library

You'll need [CMake](https://cmake.org/) and a C/C++ compiler such as
GCC, Clang or MSVC. Clone the libtopotoolbox repository

```
> git clone https://github.com/TopoToolbox/libtopotoolbox
```

and create the project's build system by entering the repository root
directory, and running cmake

```
> cd libtopotoolbox
> cmake -B build
```

Then build the project by running

```
> cmake --build build
```

By default, CMake will build the Debug version of the library, which
includes extra information for debugging purposes, but which turns off
optimizations. The Release build, which is significantly faster, can
also be built by specifying `-DCMAKE_BUILD_TYPE=Release` during the
generation step:

```
> cmake -B -DCMAKE_BUILD_TYPE=Release
```

libtopotoolbox includes a test suite that can be built alongside the
library. Turn on the `TT_BUILD_TESTS` option to build the tests as well:

```
> cmake -B build -DTT_BUILD_TESTS=ON
> cmake --build build
```

Tests can then be run:

```
> ctest --test-dir build
```

If you make any changes to the existing source files or tests, you can
run the two previous commands to verify that the library builds and
the tests pass.

### Formatting and style

Unless otherwise stated, we follow Google's [C++ style
guide](https://google.github.io/styleguide/cppguide.html).

Code formatting is enforced by `clang-format` during the automated
checks of each pull request. A `.clang-format` file is provided in the
root of the repository so that you can run `clang-format` on the
command line or set up your editor or IDE to autoformat your code.

Library code (in `src/`) should target C99 without compiler-specific
extensions. The build system supports C++11, but C++ code should only
be contributed to the library with a clear justification. Test code
(in `test/`) can be written in C or C++ as needed.

All public-facing functions exported by libtopotoolbox should be
declared in `include/topotoolbox.h` with the `TOPOTOOLBOX_API`
specifier. For example:

``` C
TOPOTOOLBOX_API int has_topotoolbox(void);
```
