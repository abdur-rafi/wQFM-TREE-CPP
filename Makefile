# Prototype build; CMake replaces this in Phase 1.
# Needs a C++20 compiler with <barrier> (Apple clang 16, GCC >= 11, clang >= 14), a JDK and Python 3.

CXXFLAGS ?= -std=c++20 -O3 -ffp-contract=off -Wall -Wextra -pthread
ifeq ($(shell uname -m),arm64)
  NATIVE ?= -mcpu=native
else
  NATIVE ?= -march=native
endif
WQFM_TREE_DIR ?= ..

all: build/gainsweep

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
