# EXECUTABLE PROJECT - bench

# === Create Executable ===
add_executable(bench "bench.cpp" "Timer.h")

# === Libraries ===
find_package(PkgConfig)
pkg_check_modules(mpfr REQUIRED IMPORTED_TARGET mpfr)
target_link_libraries(bench PRIVATE PkgConfig::mpfr)

target_link_libraries(bench PUBLIC ReferenceLambertW)
target_include_directories(bench PUBLIC "../include/")

find_package(flint REQUIRED)
target_link_libraries(bench PRIVATE "${flint_DIR}/../../flint.lib")
target_include_directories(bench PRIVATE "${flint_DIR}/../../../include")

# === Feature Enables ===
if (REFERENCEW_MSVC_STATIC_RUNTIME)
    set_property(TARGET bench PROPERTY MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
endif()

target_compile_features(bench PUBLIC cxx_std_20)