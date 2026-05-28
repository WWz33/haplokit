cmake_minimum_required(VERSION 3.10)
project(gffsub VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(gffsub
    src/gffsub.cpp
    src/selector_filter.cpp
    src/qc_validate.cpp
    src/qc_parse.cpp
    src/gtf_parser.cpp
    src/query_summary.cpp
    src/attributes.cpp
    src/annotation_index.cpp
    src/gff3_parser.cpp
    src/region.cpp
    src/annotation_filter.cpp
    src/isoform_filter.cpp
    src/annotation_output.cpp
)

target_include_directories(gffsub PRIVATE
    src
)

target_compile_options(gffsub PRIVATE
    -Wall -O2
)

include(CTest)

if(BUILD_TESTING)
    add_executable(annotation_index_smoke
        tests/annotation_index_smoke.cpp
        src/attributes.cpp
        src/annotation_index.cpp
        src/gff3_parser.cpp
        src/gtf_parser.cpp
        src/region.cpp
    )

    target_include_directories(annotation_index_smoke PRIVATE
        src
    )

    target_compile_options(annotation_index_smoke PRIVATE
        -Wall -O2
    )

    add_test(NAME annotation_index_smoke COMMAND annotation_index_smoke)

    add_executable(cli_output_attrs_smoke
        tests/cli_output_attrs_smoke.cpp
    )

    target_compile_options(cli_output_attrs_smoke PRIVATE
        -Wall -O2
    )

    add_test(NAME cli_output_attrs_smoke COMMAND cli_output_attrs_smoke $<TARGET_FILE:gffsub>)

    add_executable(cli_selector_smoke
        tests/cli_selector_smoke.cpp
    )

    target_compile_options(cli_selector_smoke PRIVATE
        -Wall -O2
    )

    add_test(NAME cli_selector_smoke COMMAND cli_selector_smoke $<TARGET_FILE:gffsub>)
endif()

# Installation
install(TARGETS gffsub DESTINATION bin)
