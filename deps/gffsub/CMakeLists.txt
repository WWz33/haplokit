cmake_minimum_required(VERSION 3.10)
project(gffsub VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find htslib
find_package(PkgConfig REQUIRED)
pkg_check_modules(HTSLIB REQUIRED htslib)

add_executable(gffsub
    src/gffsub.cpp
    src/gff3_parser.cpp
    src/gff3_filter.cpp
    src/gff3_output.cpp
)

target_include_directories(gffsub PRIVATE
    ${HTSLIB_INCLUDE_DIRS}
    src
)

target_link_libraries(gffsub PRIVATE
    ${HTSLIB_LIBRARIES}
)

target_compile_options(gffsub PRIVATE
    -Wall -O2
)

# Installation
install(TARGETS gffsub DESTINATION bin)

# Print status
message(STATUS "HTSLIB include: ${HTSLIB_INCLUDE_DIRS}")
message(STATUS "HTSLIB libraries: ${HTSLIB_LIBRARIES}")