function(haplokit_configure_htslib)
    option(HAPTOOLS_USE_VENDORED_HTSLIB "Build and link vendored htslib from deps/htslib" ON)
    set(HAPTOOLS_VENDORED_HTSLIB_DIR "${CMAKE_SOURCE_DIR}/deps/htslib" CACHE PATH "Path to vendored htslib source tree")

    if(HAPTOOLS_USE_VENDORED_HTSLIB)
        if(NOT EXISTS "${HAPTOOLS_VENDORED_HTSLIB_DIR}/Makefile")
            message(FATAL_ERROR "Vendored htslib not found at ${HAPTOOLS_VENDORED_HTSLIB_DIR}")
        endif()

        if(NOT EXISTS "${HAPTOOLS_VENDORED_HTSLIB_DIR}/htscodecs/htscodecs/htscodecs.c")
            message(FATAL_ERROR
                "Vendored htscodecs sources are missing. Populate deps/htslib/htscodecs before building."
            )
        endif()

        add_custom_command(
            OUTPUT "${HAPTOOLS_VENDORED_HTSLIB_DIR}/libhts.a"
            COMMAND "${CMAKE_MAKE_PROGRAM}" -C "${HAPTOOLS_VENDORED_HTSLIB_DIR}" lib-static -j1
            WORKING_DIRECTORY "${HAPTOOLS_VENDORED_HTSLIB_DIR}"
            COMMENT "Building vendored htslib"
            VERBATIM
        )

        add_custom_target(haplokit_htslib_build DEPENDS "${HAPTOOLS_VENDORED_HTSLIB_DIR}/libhts.a")

        add_library(haplokit_htslib STATIC IMPORTED GLOBAL)
        set_target_properties(
            haplokit_htslib
            PROPERTIES
                IMPORTED_LOCATION "${HAPTOOLS_VENDORED_HTSLIB_DIR}/libhts.a"
                INTERFACE_INCLUDE_DIRECTORIES "${HAPTOOLS_VENDORED_HTSLIB_DIR}"
        )
        add_dependencies(haplokit_htslib haplokit_htslib_build)

        set(_haplokit_htslib_target haplokit_htslib)
        set(_haplokit_htslib_include_dirs "${HAPTOOLS_VENDORED_HTSLIB_DIR}")
        set(_haplokit_htslib_extra_libs z m bz2 lzma curl pthread)
    else()
        find_package(PkgConfig REQUIRED)
        pkg_check_modules(HTSLIB REQUIRED htslib)

        add_library(haplokit_htslib INTERFACE)
        target_include_directories(haplokit_htslib INTERFACE ${HTSLIB_INCLUDE_DIRS})
        target_link_directories(haplokit_htslib INTERFACE ${HTSLIB_LIBRARY_DIRS})
        target_link_libraries(haplokit_htslib INTERFACE ${HTSLIB_LIBRARIES})
        target_compile_options(haplokit_htslib INTERFACE ${HTSLIB_CFLAGS_OTHER})

        set(_haplokit_htslib_target haplokit_htslib)
        set(_haplokit_htslib_include_dirs ${HTSLIB_INCLUDE_DIRS})
        set(_haplokit_htslib_extra_libs "")
    endif()

    set(HAPTOOLS_HTSLIB_TARGET "${_haplokit_htslib_target}" PARENT_SCOPE)
    set(HAPTOOLS_HTSLIB_INCLUDE_DIRS "${_haplokit_htslib_include_dirs}" PARENT_SCOPE)
    set(HAPTOOLS_HTSLIB_EXTRA_LIBS "${_haplokit_htslib_extra_libs}" PARENT_SCOPE)
endfunction()

