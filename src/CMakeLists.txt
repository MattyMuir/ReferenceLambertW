add_executable(ReferenceLambertW "main.cpp" "ReferenceLambertW.h" "ReferenceLambertW.cpp" "ReferenceW.h" "ReferenceW.cpp" "Interval.h")

# Link against MPFR
find_package(PkgConfig)
pkg_check_modules(mpfr REQUIRED IMPORTED_TARGET mpfr)
target_link_libraries(ReferenceLambertW PRIVATE PkgConfig::mpfr)

# Use C++20
set_property(TARGET ReferenceLambertW PROPERTY CXX_STANDARD 20)