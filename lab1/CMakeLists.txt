cmake_minimum_required(VERSION 3.16)
project(macierzowe-lab1 LANGUAGES CXX)

include(CTest)
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ==========================
# BIBLIOTEKA Z LOGIKĄ
# ==========================
add_library(multiply_binet src/multiply_binet.cpp)
target_include_directories(multiply_binet PUBLIC include)
target_compile_options(multiply_binet PRIVATE -Wall -Wextra -pedantic -O3)

add_library(multiply_strassen src/multiply_strassen.cpp)
target_include_directories(multiply_strassen PUBLIC include)
target_compile_options(multiply_strassen PRIVATE -Wall -Wextra -pedantic -O3)

# ==========================
# GŁÓWNY PROGRAM
# ==========================
add_executable(${PROJECT_NAME} src/main.cpp)
target_link_libraries(${PROJECT_NAME} PRIVATE multiply_binet multiply_strassen)

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
