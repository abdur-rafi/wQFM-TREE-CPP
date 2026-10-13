# Prototype build; CMake replaces this in Phase 1.
# Needs a C++20 compiler with <barrier> (Apple clang 16, GCC >= 11, clang >= 14), a JDK and Python 3.

CXXFLAGS ?= -std=c++20 -O3 -ffp-contract=off -Wall -Wextra -pthread
ifeq ($(shell uname -m),arm64)
  NATIVE ?= -mcpu=native
else
  NATIVE ?= -march=native
endif
WQFM_TREE_DIR ?= ..

LIB_SRC := src/newick.cpp src/forest.cpp src/scorer.cpp src/consensus.cpp src/trace.cpp src/wqfm.cpp
HEADERS := $(wildcard src/*.hpp)

all: build/wqfm-tree build/brute_check build/gainsweep

build/wqfm-tree: $(LIB_SRC) src/main.cpp $(HEADERS)
	@mkdir -p build
	$(CXX) $(CXXFLAGS) $(NATIVE) -Isrc -o $@ $(LIB_SRC) src/main.cpp

build/brute_check: $(LIB_SRC) tests/brute_check.cpp $(HEADERS)
	@mkdir -p build
	$(CXX) $(CXXFLAGS) $(NATIVE) -Isrc -o $@ $(LIB_SRC) tests/brute_check.cpp

build/gainsweep: prototype/gainsweep.cpp
	@mkdir -p build
	$(CXX) $(CXXFLAGS) $(NATIVE) -o $@ $<

java:
	tools/build_java.sh $(WQFM_TREE_DIR)

check: all java
	tools/check.sh

clean:
	rm -rf build

.PHONY: all java check clean
