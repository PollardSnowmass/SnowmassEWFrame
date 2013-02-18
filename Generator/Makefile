run: pythia_run.cxx
	g++ pythia_run.cxx -o run -I$(PYTHIA8)/include -I$(FASTJET)/include \
		-I$(LHAPDF)/include -L$(PYTHIA8)/lib -L$(FASTJET)/lib -L$(LHAPDF)/lib \
		-lpythia8 -lfastjet -lLHAPDF `root-config --libs --cflags`

clean:
	rm *.o *.so run
