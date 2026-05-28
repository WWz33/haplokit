CXX ?= g++
CXXFLAGS = -O2 -Wall -std=c++17 -pthread
PREFIX ?= /usr/local

TARGET = gffsub
SRCS = src/gffsub.cpp \
       src/selector_filter.cpp \
       src/qc_validate.cpp \
       src/qc_parse.cpp \
       src/gtf_parser.cpp \
       src/query_summary.cpp \
       src/attributes.cpp \
       src/annotation_index.cpp \
       src/gff3_parser.cpp \
       src/region.cpp \
       src/annotation_filter.cpp \
       src/isoform_filter.cpp \
       src/annotation_output.cpp
HDRS = src/gff3.hpp \
       src/annotation.hpp \
       src/qc_validate.hpp \
       src/qc_parse.hpp \
       src/gtf_parser.hpp \
       src/query_summary.hpp \
       src/selector_filter.hpp
ANNOTATION_INDEX_SMOKE = annotation_index_smoke
CLI_OUTPUT_ATTRS_SMOKE = cli_output_attrs_smoke
CLI_SELECTOR_SMOKE = cli_selector_smoke

.PHONY: all clean test install uninstall

all: $(TARGET)

$(TARGET): $(SRCS) $(HDRS)
	$(CXX) $(CXXFLAGS) -o $@ $(SRCS)

$(ANNOTATION_INDEX_SMOKE): tests/annotation_index_smoke.cpp src/attributes.cpp src/annotation_index.cpp src/gff3_parser.cpp src/gtf_parser.cpp src/region.cpp $(HDRS)
	$(CXX) $(CXXFLAGS) -Isrc -o $@ tests/annotation_index_smoke.cpp src/attributes.cpp src/annotation_index.cpp src/gff3_parser.cpp src/gtf_parser.cpp src/region.cpp

$(CLI_OUTPUT_ATTRS_SMOKE): tests/cli_output_attrs_smoke.cpp
	$(CXX) $(CXXFLAGS) -o $@ tests/cli_output_attrs_smoke.cpp

$(CLI_SELECTOR_SMOKE): tests/cli_selector_smoke.cpp
	$(CXX) $(CXXFLAGS) -o $@ tests/cli_selector_smoke.cpp

clean:
	rm -f $(TARGET) $(ANNOTATION_INDEX_SMOKE) $(CLI_OUTPUT_ATTRS_SMOKE) $(CLI_SELECTOR_SMOKE) src/*.o

test: $(TARGET) $(ANNOTATION_INDEX_SMOKE) $(CLI_OUTPUT_ATTRS_SMOKE) $(CLI_SELECTOR_SMOKE)
	./$(ANNOTATION_INDEX_SMOKE)
	./$(CLI_OUTPUT_ATTRS_SMOKE) ./$(TARGET)
	./$(CLI_SELECTOR_SMOKE) ./$(TARGET)

install: $(TARGET)
	install -d $(PREFIX)/bin
	install -m 755 $(TARGET) $(PREFIX)/bin/

uninstall:
	rm -f $(PREFIX)/bin/$(TARGET)
