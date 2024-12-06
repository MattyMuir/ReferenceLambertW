# LIBRARY PROJECT - ReferenceLambertW

# == Includes ==
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

# === Create Library ===
add_library(ReferenceLambertW "Interval.h" "ReferenceW.cpp"  "ReferenceW.h"  "halley.h" "ReferenceWf.h" "ReferenceWf.cpp" "rndutil.h" "rndutil.cpp" "Sign.h" )

# === Libraries ===
find_package(PkgConfig)
pkg_check_modules(mpfr REQUIRED IMPORTED_TARGET mpfr)
target_link_libraries(ReferenceLambertW PUBLIC PkgConfig::mpfr)

find_package(sleef REQUIRED)
target_link_libraries(ReferenceLambertW PUBLIC sleef::sleef)
target_include_directories(ReferenceLambertW PRIVATE ${SLEEF_INCLUDE_DIR})

find_package(flint REQUIRED)
target_link_libraries(ReferenceLambertW PUBLIC flint::flint)

# === Feature Enables ===
if (REFERENCEW_MSVC_STATIC_RUNTIME)
    set_property(TARGET ReferenceLambertW PROPERTY MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
endif()
target_compile_features(ReferenceLambertW PUBLIC cxx_std_20)
enable_ipo(ReferenceLambertW)
set_arch(ReferenceLambertW)
enable_strict_math(ReferenceLambertW)

# === Install ===
# Install binaries and headers
install(TARGETS ${PROJECT_NAME}
    EXPORT ${PROJECT_NAME}Config
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} # static library files
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} # includes
)

# Export includes
install(DIRECTORY "../include"
        DESTINATION ${CMAKE_INSTALL_PREFIX})
install(DIRECTORY "../src"
        DESTINATION ${CMAKE_INSTALL_PREFIX}
        FILES_MATCHING
        PATTERN "*.h")

# Configure the main config.cmake file
configure_package_config_file(
    ${CMAKE_SOURCE_DIR}/cmake/${PROJECT_NAME}Config.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}
)

# Install the export set and configuration files
install(EXPORT ${PROJECT_NAME}Config
    FILE ${PROJECT_NAME}Targets.cmake
    NAMESPACE ${PROJECT_NAME}::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}
)

install(
    FILES
        ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}
)