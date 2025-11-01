cmake_minimum_required(VERSION 3.16)
project(macierzowe-lab1 LANGUAGES CXX)

include(CTest)
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ==========================
# BIBLIOTEKA Z LOGIKĄ
# ==========================
add_library(matrix_mult
    src/multiply_binet.cpp
)
target_include_directories(matrix_mult PUBLIC include)
target_compile_options(matrix_mult PRIVATE -Wall -Wextra -pedantic)

# ==========================
# GŁÓWNY PROGRAM
# ==========================
add_executable(${PROJECT_NAME} src/main.cpp)
target_link_libraries(${PROJECT_NAME} PRIVATE matrix_mult)


# ==========================
# GOOGLE TEST
# ==========================
include(FetchContent)
FetchContent_Declare(
    googletest
    URL https://github.com/google/googletest/archive/refs/tags/v1.15.0.zip
)
FetchContent_MakeAvailable(googletest)

enable_testing()
add_subdirectory(test)
