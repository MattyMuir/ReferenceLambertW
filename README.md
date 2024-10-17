# Reference Lambert W Implementation

A robust reference implementation of the Lambert W function using MPFR.

## Building on Windows
This library depends on MPFR. On windows I recommend using `vcpkg` to manage this.

First, use `vcpkg` to install `mpfr` and `pkgconf`, and then specify the following CMake variables before configuring
- `CMAKE_TOOLCHAIN_FILE` = \{VCPKG_ROOT\}/scripts/buildsystems/vcpkg.cmake
- `PKG_CONFIG_EXECUTABLE` = \{VCPKG_ROOT\}/installed/x64-windows/tools/pkgconf/pkgconf.exe