cmake_minimum_required(VERSION 3.16)
project(macierzowe-lab1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


add_executable(${PROJECT_NAME}
    src/main.cpp
)

target_compile_options(${PROJECT_NAME} PRIVATE -Wall -Wextra -pedantic)
