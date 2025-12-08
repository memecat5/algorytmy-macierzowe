cmake_minimum_required(VERSION 3.16)
project(macierzowe-lab3)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# dla clangd
set(CMAKE_EXPORT_COMPILE_COMMANDS 1)


include(CTest)
include(FetchContent)
FetchContent_Declare(
  Eigen3
  URL https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.zip
)
FetchContent_MakeAvailable(Eigen3)

add_executable(${PROJECT_NAME} src/main.cpp)

target_link_libraries(macierzowe-lab3 PRIVATE
    stdc++exp # dla std::print w minGW, w sumie można później wywalić
    Eigen3::Eigen
    )
