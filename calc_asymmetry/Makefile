
all: testcode.cc
	g++ -O3 -c testcode.cc -I $(shell lhapdf-config --incdir) -I $(ROOTSYS)/include
	g++ -O3 testcode.o -o testcode \
  -I $(shell lhapdf-config --incdir) \
  -L $(shell lhapdf-config --libdir) -lLHAPDF \
  -Wl,--no-as-needed $(shell root-config --glibs)

clean:
	rm -rf *.o testcode *~
