CXX = g++ -std=c++11
CXXFLAGS = -c
BOOSTFLAGS = -lboost_system -lboost_filesystem

OBJDIR = ./objs/
OBJECTS = $(OBJDIR)Graph.o $(OBJDIR)hashGraph.o $(OBJDIR)hash2DGraph.o $(OBJDIR)indexGraph.o
MCOBJECTS = $(OBJDIR)Graph.o $(OBJDIR)hashGraph.o $(OBJDIR)hash2DGraph.o $(OBJDIR)indexGraph.o $(OBJDIR)indexMCGraph.o

all: testTri2D testCNGlobal testCNGlobal2D testCNLocal testTri testTrit testMC testIndexTri hashMinMax genPairs genNeoGraph

testIndexTri :$(OBJECTS) $(OBJDIR)testIndexTri.o
	$(CXX) -o testIndexTri $(OBJECTS) $(OBJDIR)testIndexTri.o $(BOOSTFLAGS)

testMC :$(MCOBJECTS) $(OBJDIR)testMC.o
	$(CXX) -o testMC $(MCOBJECTS) $(OBJDIR)testMC.o $(BOOSTFLAGS)

testTri :$(OBJECTS) $(OBJDIR)testTri.o
	$(CXX) -o testTri $(OBJECTS) $(OBJDIR)testTri.o $(BOOSTFLAGS)

testTriBatch :$(OBJECTS) $(OBJDIR)testTriBatch.o
	$(CXX) -o testTriBatch $(OBJECTS) $(OBJDIR)testTriBatch.o $(BOOSTFLAGS)

genNeoGraph :$(OBJECTS) $(OBJDIR)genNeoGraph.o
	$(CXX) -Ofast -o genNeoGraph $(OBJECTS) $(OBJDIR)genNeoGraph.o $(BOOSTFLAGS)

testTri2D :$(OBJECTS) $(OBJDIR)testTri2D.o
	$(CXX) -o testTri2D $(OBJECTS) $(OBJDIR)testTri2D.o $(BOOSTFLAGS)

genPairs :$(OBJECTS) $(OBJDIR)genPairs.o
	$(CXX) -o genPairs $(OBJECTS) $(OBJDIR)genPairs.o $(BOOSTFLAGS)

hashMinMax :$(OBJECTS) $(OBJDIR)hashMinMax.o
	$(CXX) -o hashMinMax $(OBJECTS) $(OBJDIR)hashMinMax.o $(BOOSTFLAGS)

testTrit :$(OBJECTS) $(OBJDIR)testTrit.o
	$(CXX) -o testTrit $(OBJECTS) $(OBJDIR)testTrit.o $(BOOSTFLAGS)

testCNLocal :$(OBJECTS) $(OBJDIR)testCNLocal.o
	$(CXX) -o testCNLocal $(OBJECTS) $(OBJDIR)testCNLocal.o $(BOOSTFLAGS)

testCNGlobal :$(OBJECTS) $(OBJDIR)testCNGlobal.o
	$(CXX) -o testCNGlobal $(OBJECTS) $(OBJDIR)testCNGlobal.o $(BOOSTFLAGS)

testCNGlobal2D :$(OBJECTS) $(OBJDIR)testCNGlobal2D.o
	$(CXX) -o testCNGlobal2D $(OBJECTS) $(OBJDIR)testCNGlobal2D.o $(BOOSTFLAGS)

$(OBJDIR)testMC.o : testMC.cpp indexGraph.h indexMCGraph.h hashGraph.h Graph.h utils.h
	$(CXX) $(CXXFLAGS) -o $(OBJDIR)testMC.o testMC.cpp

$(OBJDIR)testTri.o : testTri.cpp indexGraph.h hashGraph.h Graph.h utils.h
	$(CXX) $(CXXFLAGS) -o $(OBJDIR)testTri.o testTri.cpp

$(OBJDIR)testTriBatch.o : testTriBatch.cpp indexGraph.h hashGraph.h Graph.h utils.h
	$(CXX) $(CXXFLAGS) -o $(OBJDIR)testTriBatch.o testTriBatch.cpp

$(OBJDIR)genNeoGraph.o : genNeoGraph.cpp indexGraph.h hashGraph.h Graph.h utils.h
	$(CXX) $(CXXFLAGS) -o $(OBJDIR)genNeoGraph.o genNeoGraph.cpp
	
$(OBJDIR)testTri2D.o : testTri2D.cpp indexGraph.h hash2DGraph.h hashGraph.h Graph.h utils.h
	$(CXX) $(CXXFLAGS) -o $(OBJDIR)testTri2D.o testTri2D.cpp

$(OBJDIR)genPairs.o : genPairs.cpp indexGraph.h hashGraph.h Graph.h utils.h
	$(CXX) $(CXXFLAGS) -o $(OBJDIR)genPairs.o genPairs.cpp

$(OBJDIR)hashMinMax.o : hashMinMax.cpp indexGraph.h hashGraph.h Graph.h utils.h
	$(CXX) $(CXXFLAGS) -o $(OBJDIR)hashMinMax.o hashMinMax.cpp

$(OBJDIR)testTrit.o : testTrit.cpp indexGraph.h hashGraph.h Graph.h utils.h
	$(CXX) $(CXXFLAGS) -o $(OBJDIR)testTrit.o testTrit.cpp

$(OBJDIR)testIndexTri.o : testIndexTri.cpp indexGraph.h hashGraph.h Graph.h utils.h
	$(CXX) $(CXXFLAGS) -o $(OBJDIR)testIndexTri.o testIndexTri.cpp

$(OBJDIR)testCNLocal.o : testCNLocal.cpp indexGraph.h hashGraph.h Graph.h utils.h
	$(CXX) $(CXXFLAGS) -o $(OBJDIR)testCNLocal.o testCNLocal.cpp

$(OBJDIR)testCNGlobal.o : testCNGlobal.cpp indexGraph.h hashGraph.h hash2DGraph.h Graph.h utils.h
	$(CXX) $(CXXFLAGS) -o $(OBJDIR)testCNGlobal.o testCNGlobal.cpp

$(OBJDIR)testCNGlobal2D.o : testCNGlobal2D.cpp indexGraph.h hashGraph.h hash2DGraph.h Graph.h utils.h
	$(CXX) $(CXXFLAGS) -o $(OBJDIR)testCNGlobal2D.o testCNGlobal2D.cpp

$(OBJDIR)Graph.o : Graph.cpp Graph.h utils.h
	$(CXX) $(CXXFLAGS) -o $(OBJDIR)Graph.o Graph.cpp

$(OBJDIR)hashGraph.o : hashGraph.cpp hashGraph.h Graph.h utils.h
	$(CXX) $(CXXFLAGS) -o $(OBJDIR)hashGraph.o hashGraph.cpp

$(OBJDIR)hash2DGraph.o : hash2DGraph.cpp hash2DGraph.h Graph.h utils.h
	$(CXX) $(CXXFLAGS) -o $(OBJDIR)hash2DGraph.o hash2DGraph.cpp

$(OBJDIR)indexGraph.o : indexGraph.cpp indexGraph.h hashGraph.h Graph.h utils.h
	$(CXX) $(CXXFLAGS) -o $(OBJDIR)indexGraph.o indexGraph.cpp

$(OBJDIR)indexMCGraph.o : indexMCGraph.cpp indexMCGraph.h utils.h
	$(CXX) $(CXXFLAGS) -o $(OBJDIR)indexMCGraph.o indexMCGraph.cpp

clean:
	rm -rf $(OBJDIR)*