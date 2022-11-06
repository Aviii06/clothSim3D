list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/geometry-central/cmake")
find_package (Eigen3 3.4 QUIET)
if(EIGEN3_FOUND)
    message(STATUS "--    success. find_package() succeeded with Eigen at ${EIGEN3_INCLUDE_DIR}")
# The find_package() failed, manually download the package
else()
    message(STATUS "--    Could not find Eigen anywhere. Downloading it (or using cached download).")

    if(${CMAKE_VERSION} VERSION_LESS 3.11)
        include(FetchContentLocal)
    else()
        include(FetchContent)
    endif()
    FetchContent_Declare(
        eigen
        GIT_REPOSITORY  https://gitlab.com/libeigen/eigen.git
        GIT_TAG         3.4.0
        SOURCE_DIR      "${CMAKE_CURRENT_BINARY_DIR}/eigen-src"
        BINARY_DIR      "${CMAKE_CURRENT_BINARY_DIR}/eigen-build"
    )
    FetchContent_GetProperties(eigen)
    if(NOT eigen_POPULATED)
        FetchContent_Populate(eigen)
        include(EigenChecker)
        eigen3checker(${eigen_SOURCE_DIR} 3.4)
    endif()
    message(STATUS "--    success. Downloaded Eigen in ${eigen_SOURCE_DIR}")
endif()

add_library (eigen INTERFACE)
target_include_directories(eigen INTERFACE
    $<BUILD_INTERFACE:${EIGEN3_INCLUDE_DIR}>
    $<INSTALL_INTERFACE:${INCLUDE_INSTALL_DIR}>
)

add_subdirectory(geometry-central)

set(BUILD_SHARED_LIBS OFF)
unset(BUILD_SHARED_LIBS)
