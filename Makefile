# Build haplokit's C++ backends directly with make.
# Requires the deps/gffsub and deps/htslib submodules:
#   git submodule update --init --recursive

CXX       ?= c++
PYTHON    ?= python3
BUILD_DIR ?= build-wsl

CXXFLAGS  ?= -std=c++17 -O2 -DNDEBUG -Wall
INCLUDES  := -Isrc/cpp -Ideps/gffsub/src -Ideps/nlohmann -Ideps/htslib
# Mirrors htslib's non-configure default link set, plus pthreads.
LDLIBS    := -lz -lm -lbz2 -llzma -lcurl -lpthread

HTSLIB_LIB          := deps/htslib/libhts.a
HTSCODECS_VERSION_H := deps/htslib/htscodecs/htscodecs/version.h

# The network backend compiles an __AVX2__ Hamming path when the compiler
# accepts -mavx2; a runtime CPU check picks the implementation.
AVX2_FLAG := $(shell $(CXX) -mavx2 -x c++ -fsyntax-only /dev/null 2>/dev/null && echo -mavx2)

CORE_SRCS := \
    src/cpp/selector.cpp \
    src/cpp/vcf_reader.cpp \
    src/cpp/view_backend.cpp \
    src/cpp/gff_annotator.cpp \
    deps/gffsub/src/feature_types.cpp \
    deps/gffsub/src/string_utils.cpp \
    deps/gffsub/src/record.cpp \
    deps/gffsub/src/gff3_parser.cpp \
    deps/gffsub/src/gtf_parser.cpp \
    deps/gffsub/src/attributes.cpp \
    deps/gffsub/src/annotation_index.cpp \
    deps/gffsub/src/region.cpp \
    deps/gffsub/src/annotation_filter.cpp \
    deps/gffsub/src/isoform_filter.cpp \
    deps/gffsub/src/annotation_output.cpp

NETWORK_SRCS := \
    src/cpp/network/graph.cpp \
    src/cpp/network/hamming_simd.cpp \
    src/cpp/network/tcs_algorithm.cpp \
    src/cpp/network/msn_algorithm.cpp \
    src/cpp/network/mjn_algorithm.cpp \
    src/cpp/network/network_backend.cpp

# Object names drop the .cpp suffix (src/cpp/main.cpp becomes
# build-wsl/obj/src/cpp/main.o) so the %.o: %.cpp pattern rule stems work.
CORE_OBJS    := $(addprefix $(BUILD_DIR)/obj/,$(addsuffix .o,$(basename $(CORE_SRCS))))
NETWORK_OBJS := $(addprefix $(BUILD_DIR)/obj/,$(addsuffix .o,$(basename $(NETWORK_SRCS))))

BINARIES := \
    $(BUILD_DIR)/haplokit_cpp \
    $(BUILD_DIR)/haplokit_network_backend \
    $(BUILD_DIR)/vcf_reader_test \
    $(BUILD_DIR)/view_backend_test \
    $(BUILD_DIR)/network_test

.PHONY: all build test test-cpp test-python check-submodules clean

all: build

build: $(BINARIES)

$(BINARIES): | check-submodules

check-submodules:
	@test -f deps/htslib/Makefile || { echo "deps/htslib submodule missing. Run: git submodule update --init --recursive" >&2; exit 1; }
	@test -f deps/htslib/htscodecs/htscodecs/htscodecs.c || { echo "htscodecs sources missing. Run: git submodule update --init --recursive" >&2; exit 1; }
	@test -f deps/gffsub/src/gff3.hpp || { echo "deps/gffsub submodule missing. Run: git submodule update --init --recursive" >&2; exit 1; }

$(BUILD_DIR)/obj/%.o: %.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

# Network sources include their headers as "graph.h" without a directory
# prefix (tests/cpp/network_test.cpp does the same).
NETWORK_INCLUDES := $(INCLUDES) -Isrc/cpp/network
$(BUILD_DIR)/obj/src/cpp/network/%.o: %.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(NETWORK_INCLUDES) $(AVX2_FLAG) -c $< -o $@

$(BUILD_DIR)/obj/tests/cpp/network_test.o: tests/cpp/network_test.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(NETWORK_INCLUDES) $(AVX2_FLAG) -c $< -o $@

# htslib's Makefile rewrites htscodecs/htscodecs/version.h from `git describe`
# whenever htscodecs/.git exists, and that rule has a force prerequisite, so it
# runs on every build. A shallow submodule clone has no tags and the version
# becomes "1.6.6-g<hash>" instead of the configured one. Stamp the
# configure.ac version before and after the build so the result does not
# depend on clone depth.
$(HTSLIB_LIB):
	@vers=`sed -n 's/^AC_INIT(htscodecs, *\([0-9.]*\).*/\1/p' deps/htslib/htscodecs/configure.ac`; \
	test -n "$$vers" || { echo "cannot determine htscodecs version" >&2; exit 1; }; \
	printf '#define HTSCODECS_VERSION_TEXT "%s"\n' "$$vers" > $(HTSCODECS_VERSION_H)
	$(MAKE) -C deps/htslib lib-static -j1
	@vers=`sed -n 's/^AC_INIT(htscodecs, *\([0-9.]*\).*/\1/p' deps/htslib/htscodecs/configure.ac`; \
	printf '#define HTSCODECS_VERSION_TEXT "%s"\n' "$$vers" > $(HTSCODECS_VERSION_H)

$(BUILD_DIR)/haplokit_cpp: $(BUILD_DIR)/obj/src/cpp/main.o $(CORE_OBJS) $(HTSLIB_LIB)
	$(CXX) $(CXXFLAGS) -o $@ $(filter %.o,$^) $(HTSLIB_LIB) $(LDLIBS)

$(BUILD_DIR)/vcf_reader_test: $(BUILD_DIR)/obj/tests/cpp/vcf_reader_test.o $(CORE_OBJS) $(HTSLIB_LIB)
	$(CXX) $(CXXFLAGS) -o $@ $(filter %.o,$^) $(HTSLIB_LIB) $(LDLIBS)

$(BUILD_DIR)/view_backend_test: $(BUILD_DIR)/obj/tests/cpp/view_backend_test.o $(CORE_OBJS) $(HTSLIB_LIB)
	$(CXX) $(CXXFLAGS) -o $@ $(filter %.o,$^) $(HTSLIB_LIB) $(LDLIBS)

$(BUILD_DIR)/haplokit_network_backend: $(BUILD_DIR)/obj/src/cpp/network/main.o $(NETWORK_OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(filter %.o,$^)

$(BUILD_DIR)/network_test: $(BUILD_DIR)/obj/tests/cpp/network_test.o $(NETWORK_OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(filter %.o,$^)

test: test-cpp test-python

test-cpp: build
	./$(BUILD_DIR)/vcf_reader_test
	./$(BUILD_DIR)/view_backend_test
	./$(BUILD_DIR)/network_test

test-python: build
	HAPLOKIT_CPP_BIN=$$PWD/$(BUILD_DIR)/haplokit_cpp $(PYTHON) -m pytest -q tests/python

clean:
	rm -rf $(BUILD_DIR)
	-$(MAKE) -C deps/htslib clean
