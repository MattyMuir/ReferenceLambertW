# EXECUTABLE PROJECT - example

# === Create Executable ===
add_executable(example "example.cpp" "Timer.h")

# === Libraries ===
find_package(PkgConfig)
pkg_check_modules(mpfr REQUIRED IMPORTED_TARGET mpfr)
target_link_libraries(example PRIVATE PkgConfig::mpfr)

target_link_libraries(example PUBLIC ReferenceLambertW)
target_include_directories(example PUBLIC "../include/")

# === Feature Enables ===
if (REFERENCEW_MSVC_STATIC_RUNTIME)
    set_property(TARGET example PROPERTY MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
endif()
target_compile_features(example PUBLIC cxx_std_20)