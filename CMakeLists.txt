cmake_minimum_required(VERSION 3.16)

project(field_d_star_algorithm LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(PkgConfig REQUIRED)
pkg_check_modules(SDL2 REQUIRED IMPORTED_TARGET sdl2)

add_library(field_d_star_planner
    cpp/field_d_star_planner.cpp
)

target_include_directories(field_d_star_planner
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/cpp
)

target_compile_features(field_d_star_planner PUBLIC cxx_std_17)

add_executable(field_d_star_gui
    cpp/main.cpp
)

target_link_libraries(field_d_star_gui
    PRIVATE
        field_d_star_planner
        PkgConfig::SDL2
)

target_compile_features(field_d_star_gui PRIVATE cxx_std_17)
