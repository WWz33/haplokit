function(haplokit_configure_htslib)
    option(HAPLOKIT_USE_VENDORED_HTSLIB "Build and link vendored htslib from deps/htslib" ON)
    set(HAPLOKIT_VENDORED_HTSLIB_DIR "${CMAKE_SOURCE_DIR}/deps/htslib" CACHE PATH "Path to vendored htslib source tree")

    if(HAPLOKIT_USE_VENDORED_HTSLIB)
        if(NOT EXISTS "${HAPLOKIT_VENDORED_HTSLIB_DIR}/Makefile")
            message(FATAL_ERROR
                "htslib submodule not found at ${HAPLOKIT_VENDORED_HTSLIB_DIR}. "
                "Run: git submodule update --init --recursive"
            )
        endif()

        if(NOT EXISTS "${HAPLOKIT_VENDORED_HTSLIB_DIR}/htscodecs/htscodecs/htscodecs.c")
            message(FATAL_ERROR
                "htscodecs sources are missing. Run: git submodule update --init --recursive"
            )
        endif()

        # htslib's Makefile rewrites htscodecs/htscodecs/version.h from
        # `git describe` whenever htscodecs/.git exists, and its rule carries a
        # `force` prerequisite so it runs on every build. A shallow submodule
        # clone has no tags, so `git describe` falls back to a bare hash and the
        # version becomes "1.6.6-g<hash>" instead of the configured version.
        # Stamp the configure.ac version both before the build (so the header
        # exists for a plain configure) and after it (so a rewrite by make is
        # undone), keeping the value independent of clone depth.
        set(_haplokit_htscodecs_version_h "${HAPLOKIT_VENDORED_HTSLIB_DIR}/htscodecs/htscodecs/version.h")
        file(READ "${HAPLOKIT_VENDORED_HTSLIB_DIR}/htscodecs/configure.ac" _haplokit_htscodecs_configure_ac)
        string(REGEX MATCH "AC_INIT\\(htscodecs, *([0-9]+\\.[0-9]+\\.[0-9]+)\\)" _haplokit_htscodecs_version_match "${_haplokit_htscodecs_configure_ac}")
        if(NOT CMAKE_MATCH_1)
            message(FATAL_ERROR "Unable to determine htscodecs version from configure.ac")
        endif()
        set(_haplokit_htscodecs_version_stamp "${CMAKE_BINARY_DIR}/htscodecs-version.h")
        file(WRITE "${_haplokit_htscodecs_version_stamp}" "#define HTSCODECS_VERSION_TEXT \"${CMAKE_MATCH_1}\"\n")
        file(COPY_FILE "${_haplokit_htscodecs_version_stamp}" "${_haplokit_htscodecs_version_h}" ONLY_IF_DIFFERENT)

        add_custom_command(
            OUTPUT "${HAPLOKIT_VENDORED_HTSLIB_DIR}/libhts.a"
            COMMAND "${CMAKE_MAKE_PROGRAM}" -C "${HAPLOKIT_VENDORED_HTSLIB_DIR}" lib-static -j1
            COMMAND "${CMAKE_COMMAND}" -E copy_if_different
                "${_haplokit_htscodecs_version_stamp}" "${_haplokit_htscodecs_version_h}"
            WORKING_DIRECTORY "${HAPLOKIT_VENDORED_HTSLIB_DIR}"
            COMMENT "Building vendored htslib"
            VERBATIM
        )

        add_custom_target(haplokit_htslib_build DEPENDS "${HAPLOKIT_VENDORED_HTSLIB_DIR}/libhts.a")

        add_library(haplokit_htslib STATIC IMPORTED GLOBAL)
        set_target_properties(
            haplokit_htslib
            PROPERTIES
                IMPORTED_LOCATION "${HAPLOKIT_VENDORED_HTSLIB_DIR}/libhts.a"
                INTERFACE_INCLUDE_DIRECTORIES "${HAPLOKIT_VENDORED_HTSLIB_DIR}"
        )
        add_dependencies(haplokit_htslib haplokit_htslib_build)

        set(_haplokit_htslib_target haplokit_htslib)
        set(_haplokit_htslib_include_dirs "${HAPLOKIT_VENDORED_HTSLIB_DIR}")
        find_package(Threads REQUIRED)

        set(_haplokit_htslib_extra_libs "")
        foreach(_haplokit_lib IN ITEMS z m bz2 lzma curl)
            string(TOUPPER "${_haplokit_lib}" _haplokit_lib_upper)
            find_library(HAPLOKIT_${_haplokit_lib_upper}_LIBRARY NAMES "${_haplokit_lib}")
            if(HAPLOKIT_${_haplokit_lib_upper}_LIBRARY)
                list(APPEND _haplokit_htslib_extra_libs "${HAPLOKIT_${_haplokit_lib_upper}_LIBRARY}")
            else()
                list(APPEND _haplokit_htslib_extra_libs "${_haplokit_lib}")
            endif()
        endforeach()
        list(APPEND _haplokit_htslib_extra_libs Threads::Threads)
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

    set(HAPLOKIT_HTSLIB_TARGET "${_haplokit_htslib_target}" PARENT_SCOPE)
    set(HAPLOKIT_HTSLIB_INCLUDE_DIRS "${_haplokit_htslib_include_dirs}" PARENT_SCOPE)
    set(HAPLOKIT_HTSLIB_EXTRA_LIBS "${_haplokit_htslib_extra_libs}" PARENT_SCOPE)
endfunction()
