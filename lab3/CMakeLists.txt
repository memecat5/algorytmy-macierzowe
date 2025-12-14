cmake_minimum_required(VERSION 3.14)
project(SvdCompressionSpectra)

set(CMAKE_CXX_STANDARD 17)

include(FetchContent)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# 1. Pobieranie Eigen
message(STATUS "Pobieranie Eigen...")
FetchContent_Declare(
  Eigen3
  GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
  GIT_TAG        3.4.0
)
FetchContent_MakeAvailable(Eigen3)

# 2. Pobieranie Spectra
message(STATUS "Pobieranie Spectra...")
FetchContent_Declare(
  Spectra
  GIT_REPOSITORY https://github.com/yixuan/spectra.git
  GIT_TAG        v1.0.1 # Stabilna wersja
)
FetchContent_MakeAvailable(Spectra)

add_executable(SvdCompression
    src/main.cpp
    src/compress.cpp
    src/partial_svd.cpp
    src/bitmap_utils.cpp
)

# Dodajemy sciezki naglowkowe
target_include_directories(SvdCompression PRIVATE 
    ${eigen3_SOURCE_DIR}
    ${spectra_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/include
)

# Linkowanie (chociaz to header-only, CMake lubi miec targety)
target_link_libraries(SvdCompression PRIVATE Eigen3::Eigen)

target_compile_options(SvdCompression PRIVATE -O3 -march=native)
