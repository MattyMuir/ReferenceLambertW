add_executable(ReferenceLambertW "main.cpp" "ReferenceLambertW.h" "ReferenceLambertW.cpp")

# I'm using vcpkg to link against mpfr on Windows. This requires installing mpfr and pkgconf and setting a few CMake variables
# CMAKE_TOOLCHAIN_FILE = (VCPKG_ROOT)/scripts/buildsystems/vcpkg.cmake
# PKG_CONFIG_EXECUTABLE = (VCPKG_ROOT)/installed/x64-windows/tools/pkgconf/pkgconf.exe
find_package(PkgConfig)
pkg_check_modules(mpfr REQUIRED IMPORTED_TARGET mpfr)
target_link_libraries(ReferenceLambertW PRIVATE PkgConfig::mpfr)