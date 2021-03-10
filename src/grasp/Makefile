CXX = clang++
VERSION = 0.02
CXXFLAGS = -O3 -Wall -I. -std=c++14
LDFLAGS = -lboost_iostreams -lboost_filesystem -lboost_system
TARGETS = grasp
OBJ = gspan.o  main.o

all: $(TARGETS)

$(TARGETS): $(OBJ)
	$(CXX) $(LDFLAGS) -o $(TARGETS) $(OBJ)

clean:
	rm -f $(OBJ) $(TARGETS) *~


