This library uses git submodules to manage dependencies, so after cloning you should run `git submodule update --init --recursive` to pull down dependencies. The only (**optional**) dependency which is not included as a submodule is Suitesparse, which can be installed via a package manager on Linux/MacOS, or manually with significant effort on Windows.

Geometry-central itself does not build any excutables, it is a library upon which you create other programs. For a simple example program and build system using geometry-central, see TODO.

