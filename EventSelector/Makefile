CXX=g++
CXXFLAGS=

run: main.o EventSelector.o
	$(CXX) $(CFLAGS) main.o -o run `root-config --libs` -lProof -lProofPlayer

main.o: main.cxx
	$(CXX) $(CFLAGS) main.cxx -c -o main.o `root-config --cflags`

EventSelector.o: EventSelector.cxx
	$(CXX) -fsyntax-only EventSelector.cxx -c -o EventSelector.o `root-config --cflags`
