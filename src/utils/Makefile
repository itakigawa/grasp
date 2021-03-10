CXX = clang++
VERSION = 0.01
CXXFLAGS = -O3 -g -Wall -I. -std=c++14
LDFLAGS = -L./ -lboost_iostreams 
TARGETS = evaluator finder predictor

all: $(TARGETS)

evaluator: evaluator.cpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o evaluator evaluator.cpp

finder: finder.cpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o finder finder.cpp

predictor: predictor.cpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o predictor predictor.cpp

clean:
	rm -f $(TARGETS) *~


