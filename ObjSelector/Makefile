CXX=g++
CXXFLAGS=

run: main.o ObjSelector.o
	$(CXX) $(CFLAGS) main.o -o run `root-config --libs` -lProof -lProofPlayer

main.o: main.cxx
	$(CXX) $(CFLAGS) main.cxx -c -o main.o `root-config --cflags`

ObjSelector.o: ObjSelector.cxx
	$(CXX) -fsyntax-only ObjSelector.cxx -c -o ObjSelector.o `root-config --cflags`
