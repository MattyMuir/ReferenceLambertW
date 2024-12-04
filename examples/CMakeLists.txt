# EXECUTABLE PROJECT - example

# === Create Executable ===
add_executable(example "example.cpp")
add_executable(ratapprox "ratapprox.cpp")

# === Libraries ===
find_package(PkgConfig)
pkg_check_modules(mpfr REQUIRED IMPORTED_TARGET mpfr)
target_link_libraries(example PRIVATE PkgConfig::mpfr)
target_link_libraries(ratapprox PRIVATE PkgConfig::mpfr)

target_link_libraries(example PUBLIC ReferenceLambertW)
target_link_libraries(ratapprox PUBLIC ReferenceLambertW)
target_include_directories(example PUBLIC "../include/")
target_include_directories(ratapprox PUBLIC "../include/")

find_package(flint REQUIRED)
target_link_libraries(example PRIVATE "${flint_DIR}/../../flint.lib")
target_link_libraries(ratapprox PRIVATE "${flint_DIR}/../../flint.lib")
target_include_directories(example PRIVATE "${flint_DIR}/../../../include")
target_include_directories(ratapprox PRIVATE "${flint_DIR}/../../../include")

# === Feature Enables ===
if (REFERENCEW_MSVC_STATIC_RUNTIME)
    set_property(TARGET example PROPERTY MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
    set_property(TARGET ratapprox PROPERTY MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
endif()

target_compile_features(example PUBLIC cxx_std_20)
target_compile_features(ratapprox PUBLIC cxx_std_20)