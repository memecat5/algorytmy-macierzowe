cmake_minimum_required(VERSION 3.10)

# Nazwa projektu
project(HierarchicalMatrixCompression)

# Ustawienie standardu C++ (C++17 jest bezpiecznym wyborem dla nowoczesnych bibliotek)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ------------------------------------------------------------------------------
# 1. Konfiguracja biblioteki Eigen3
# ------------------------------------------------------------------------------
# Próba znalezienia zainstalowanej biblioteki Eigen3
find_package(Eigen3 3.3 QUIET NO_MODULE)

if(NOT Eigen3_FOUND)
    message(STATUS "Eigen3 not found in system. Attempting to fetch content...")
    include(FetchContent)
    FetchContent_Declare(
        Eigen3
        GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
        GIT_TAG        3.4.0
    )
    FetchContent_MakeAvailable(Eigen3)
    # Po pobraniu przez FetchContent, Eigen jest dostępny jako target 'Eigen3::Eigen'
endif()

# ------------------------------------------------------------------------------
# 2. Definicja pliku wykonywalnego
# ------------------------------------------------------------------------------
# Zakładam, że Twój plik z kodem nazywa się main.cpp
add_executable(HMatrixCompression src/main.cpp)

# ------------------------------------------------------------------------------
# 3. Linkowanie i Include
# ------------------------------------------------------------------------------
# Linkowanie Eigen
target_link_libraries(HMatrixCompression PRIVATE Eigen3::Eigen)

# Dodanie bieżącego katalogu do include (żeby widział stb_image.h w tym samym folderze)
target_include_directories(HMatrixCompression PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

# Opcjonalnie: optymalizacje dla kompilatora (dla szybszego działania SVD)
if(MSVC)
    target_compile_options(HMatrixCompression PRIVATE /O2 /fp:fast)
else()
    target_compile_options(HMatrixCompression PRIVATE -O3 -march=native)
endif()