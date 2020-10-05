CXXFLAGS ?= -std=c++17 -g -O2 -march=native

all : test_data particles

clean:
	-rm -f test_data particles

.PHONY: all clean

test_data : test_data.cc
	g++ $(CXXFLAGS) -o $@ $<

particles : particles.cc
	mpicxx $(CXXFLAGS) -o $@ -l p4est -l sc $<
