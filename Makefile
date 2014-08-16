CXX = g++
CXXFLAGS = -std=c++11 -g -Wall -Wextra -O3

ising-3d-QA: ising-3d-QA.cc
	$(CXX) $(CXXFLAGS) -Isrc -o $@ $^

clean:
	rm -rf ising-3d-QA result*.txt

