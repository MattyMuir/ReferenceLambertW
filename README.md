# Reference Lambert W Implementation
A robust reference implementation of the Lambert W function.

## Building on Windows
This library depends on flint, MPFR and sleef. Since MPFR uses a make buildsystem, I recommend using `vcpkg` to manage this on Windows.

- First, use `vcpkg` to install `mpfr` and `pkgconf`
- Download, build, and install flint and sleef (from source or using vcpkg)
- Make sure to specify the following CMake variables before configuring
	- `CMAKE_TOOLCHAIN_FILE` = \{VCPKG_ROOT\}/scripts/buildsystems/vcpkg.cmake
	- `PKG_CONFIG_EXECUTABLE` = \{VCPKG_ROOT\}/installed/x64-windows/tools/pkgconf/pkgconf.exe