# EXECUTABLE PROJECT - tests

# === Create Executable ===
add_executable(tests "tests.cpp")

# === Libraries ===
find_package(PkgConfig)
pkg_check_modules(mpfr REQUIRED IMPORTED_TARGET mpfr)
target_link_libraries(tests PRIVATE PkgConfig::mpfr)

target_link_libraries(tests PUBLIC ReferenceLambertW)
target_include_directories(tests PUBLIC "../include/")

find_package(flint REQUIRED)
target_link_libraries(tests PRIVATE flint::flint)

# === Feature Enables ===
if (REFERENCEW_MSVC_STATIC_RUNTIME)
    set_property(TARGET tests PROPERTY MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
endif()
target_compile_features(tests PUBLIC cxx_std_20)
enable_ipo(tests)
set_arch(tests)

target_compile_features(tests PUBLIC cxx_std_20)