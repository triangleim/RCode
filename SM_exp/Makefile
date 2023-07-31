CXX = g++ -std=c++11
# EXEFLAGS = -O0
CXXFLAGS = -O0 -c
BOOSTFLAGS = -lboost_system -lboost_filesystem

SOURCES = $(shell find . | grep -e ".hpp")
OBJECTS = $(SOURCES:%.hpp=%.o)

all: getData testSM
getData: getData.cpp
	$(CXX) -O3 -o getData getData.cpp $(BOOSTFLAGS)

testSM: $(OBJECTS) testSM.cpp
	$(CXX) $^ -O0 -o $@ $(BOOSTFLAGS)

%.o: %.cpp %.hpp util.hpp
	$(CXX) $(CXXFLAGS) $< -o $@

clean:
	rm -rf getData testSM $(OBJECTS)
