CXX = g++ -O3 -std=c++11

all: globalAssign globalAssignLazy localAssignLazy localAssign localAssignOrder localAssignDegree localAssignAcc globalAssign2D

globalAssign : globalAssign.cpp
	$(CXX) -o globalAssign globalAssign.cpp
globalAssignLazy : globalAssignLazy.cpp
	$(CXX) -o globalAssignLazy globalAssignLazy.cpp
globalAssign2D : globalAssign2D.cpp
	$(CXX) -o globalAssign2D globalAssign2D.cpp
localAssign : localAssign.cpp
	$(CXX) -o localAssign localAssign.cpp
localAssignLazy : localAssignLazy.cpp
	$(CXX) -o localAssignLazy localAssignLazy.cpp
localAssignOrder : localAssignOrder.cpp
	$(CXX) -o localAssignOrder localAssignOrder.cpp
localAssignDegree : localAssignDegree.cpp
	$(CXX) -o localAssignDegree localAssignDegree.cpp
localAssignAcc : localAssignAcc.cpp
	$(CXX) -o localAssignAcc localAssignAcc.cpp
clean:
	rm -rf globalAssign globalAssignLazy localAssignLazy localAssign localAssignOrder localAssignDegree localAssignAcc globalAssign2D