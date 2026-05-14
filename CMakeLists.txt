cmake_minimum_required(VERSION 3.22)

project(haplokit LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include(CTest)
enable_testing()

include(cmake/htslib.cmake)
haplokit_configure_htslib()

add_library(haplokit_core
    src/cpp/selector.cpp
    src/cpp/vcf_reader.cpp
    src/cpp/view_backend.cpp
    src/cpp/gff_annotator.cpp
    deps/gffsub/src/gff3_parser.cpp
    deps/gffsub/src/gff3_filter.cpp
    deps/gffsub/src/gff3_output.cpp
)

target_include_directories(haplokit_core
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/src/cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/deps/gffsub/src
        ${CMAKE_CURRENT_SOURCE_DIR}/deps/nlohmann
        ${HAPTOOLS_HTSLIB_INCLUDE_DIRS}
)

target_link_libraries(haplokit_core PUBLIC ${HAPTOOLS_HTSLIB_TARGET} ${HAPTOOLS_HTSLIB_EXTRA_LIBS})
add_dependencies(haplokit_core haplokit_htslib_build)

add_executable(haplokit_cpp
    src/cpp/main.cpp
)

target_link_libraries(haplokit_cpp PRIVATE haplokit_core)

if(BUILD_TESTING)
    add_executable(vcf_reader_test
        tests/cpp/vcf_reader_test.cpp
    )

    target_link_libraries(vcf_reader_test PRIVATE haplokit_core)

    add_executable(view_backend_test
        tests/cpp/view_backend_test.cpp
    )

    target_link_libraries(view_backend_test PRIVATE haplokit_core)

    add_test(
        NAME vcf_reader_test
        COMMAND vcf_reader_test
    )

    add_test(
        NAME view_backend_test
        COMMAND view_backend_test
    )

    set_tests_properties(
        vcf_reader_test view_backend_test
        PROPERTIES
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    )
endif()
